import os
import re
import sys
import copy
import argparse
import random
from itertools import zip_longest, groupby
from dataclasses import dataclass

from pysam import AlignmentFile
from xopen import xopen

CIGAR_MATCH = 0  # M
CIGAR_INSERTION = 1  # I
CIGAR_DELETION = 2  # D
CIGAR_SOFTCLIP = 4  # S


class Scores:
    match = 2
    mismatch = 8
    gap_open = 12
    gap_extension = 1
    end_bonus = 10

@dataclass
class ReferenceInterval:
    name: str
    start: int
    end: int


def only_r1_iter(bam_iter):
    """
    Use only the R1 reads in the paired-end input BAM and present the records
    as if they were from a single-end BAM.
    """
    for record in bam_iter:
        if not record.is_read1:
            continue
        # Remove PAIRED,PROPER_PAIR,MUNMAP,MREVERSE,READ1,READ2 flags
        record.flag = record.flag & ~235
        yield record


def pick_random_primary_single_end_iter(bam_iter):
    """
    Pick one of multiple primary alignments randomly
    """
    for query_name, records in groupby(bam_iter, lambda record: record.query_name):
        records = list(records)
        yield random.choice(records)


def pick_random_primary_paired_end_iter(bam_iter):
    for query_name, records in groupby(bam_iter, lambda record: record.query_name[:-2]):
        records = list(records)
        r1_records = [record for record in records if record.query_name.endswith("/1")]
        r2_records = {(record.reference_name, record.reference_start): record for record in records if record.query_name.endswith("/2")}

        both = [
            (r1, r2_records.get((r1.next_reference_name, r1.next_reference_start)))
            for r1 in r1_records
        ]
        # TODO
        # if one of the reads is unmapped, we skip!
        both = [(r1, r2) for (r1, r2) in both if r2 is not None]
        if both:
            yield from random.choice(both)


def read_alignments(bam_path, only_r1: bool):
    bam_path = AlignmentFile(bam_path, check_sq=False)
    read_positions = {}

    bam_iter = bam_path.fetch(until_eof=True)
    if only_r1:
        bam_iter = only_r1_iter(bam_iter)
    for read in bam_iter:
        if read.is_secondary:
            continue
        if read.flag == 0 or read.flag == 16:  # single end
            query_name = read.query_name
            if not query_name.endswith("/1"):
                query_name += "/1"
            if read.reference_end is None:
                read_positions[query_name] = False
            else:
                read_positions[query_name] = ReferenceInterval(
                    read.reference_name,
                    read.reference_start,
                    read.reference_end,
                )
        elif read.is_paired:
            query_name = read.query_name
            if read.is_read1 and not query_name.endswith("/1"):
                query_name += "/1"
            if read.is_read2 and not query_name.endswith("/2"):
                query_name += "/2"

            if read.is_unmapped:
                read_positions[query_name] = False
            else:
                read_positions[query_name] = ReferenceInterval(
                    read.reference_name,
                    read.reference_start,
                    read.reference_end,
                )
        elif read.is_unmapped:  # single and unmapped
            assert not read.is_paired
            query_name = read.query_name
            if not query_name.endswith("/1"):
                query_name += "/1"
            read_positions[query_name] = False

    return read_positions


def read_paf(paf_file):
    read_positions = {}  # query_name -> ReferenceInterval
    mapped_to_multiple_pos = 0
    for line in xopen(paf_file):
        vals = line.split()
        query_name, reference_name, reference_start, reference_end = (
            vals[0],
            vals[5],
            int(vals[7]),
            int(vals[8]),
        )

        if query_name in read_positions:
            mapped_to_multiple_pos += 1
            continue
        else:
            read_positions[query_name] = ReferenceInterval(reference_name, reference_start, reference_end)
    return read_positions, mapped_to_multiple_pos


def overlap(q_a, q_b, p_a, p_b):
    assert q_a <= q_b and p_a <= p_b
    # if (q_a == q_b) or (p_a == p_b):
    #     print("Cigar bug")
    return (
        (p_a <= q_a <= p_b)
        or (p_a <= q_b <= p_b)
        or (q_a <= p_a <= q_b)
        or (q_a <= p_b <= q_b)
    )


def jaccard_overlap(a_start, a_end, b_start, b_end):
    assert a_start < a_end
    assert b_start < b_end

    intersect = min(a_end, b_end) - max(a_start, b_start)
    if intersect < 0:
        return 0
    union = max(a_end, b_end) - min(a_start, b_start)
    result = intersect / union
    assert 0 <= result <= 1.0
    return result


assert jaccard_overlap(5, 10, 5, 10) == 1
assert jaccard_overlap(10, 20, 20, 30) == 0
assert jaccard_overlap(20, 30, 10, 20) == 0
assert jaccard_overlap(0, 4, 1, 3) == 0.5
assert jaccard_overlap(1, 3, 0, 4) == 0.5


def get_stats(truth, predicted):
    nr_total = len(truth)
    unaligned = 0
    nr_aligned = 0
    over_mapped = 0
    correct = 0
    correct_jaccard = 0.0
    for query_name in predicted:
        if not truth[query_name]:
            over_mapped += 1
            continue
        if not predicted[query_name]:
            unaligned += 1
            continue

        nr_aligned += 1

        predicted_interval = predicted[query_name]
        true_interval = truth[query_name]

        if predicted_interval.name == true_interval.name:
            if overlap(
                predicted_interval.start, predicted_interval.end, true_interval.start, true_interval.end
            ):
                correct += 1
            correct_jaccard += jaccard_overlap(
                predicted_interval.start, predicted_interval.end, true_interval.start, true_interval.end
            )

    aligned_percentage = 100 * nr_aligned / nr_total
    accuracy = 100 * correct / nr_total
    return aligned_percentage, accuracy, over_mapped, 100 * correct_jaccard / nr_total


def parse_int_or_not(s):
    try:
        return int(s)
    except ValueError:
        return s


def recompute_alignment_score(segment, scores):
    md_tag = segment.get_tag("MD")
    md = [
        parse_int_or_not(e)
        for e in re.split("([A-Z]|\^[A-Z]+|[0-9]+)", md_tag)
        if e != ""
    ]
    score = 0
    # To compute the score, it is not necessary to reconstruct the full
    # alignment.
    # - The numbers in the MD tag give us the number of matches
    # - The letters in the MD tag give us the number of mismatches

    for item in md:
        if isinstance(item, int):
            # Match
            score += item * scores.match
        elif not item.startswith("^"):
            # Mismatch
            assert len(item) == 1
            score -= scores.mismatch

    # - I or D operations in the CIGAR string give us indels and their lengths
    cigartuples = segment.cigartuples
    for op, length in cigartuples:
        if op == CIGAR_INSERTION or op == CIGAR_DELETION:
            score -= scores.gap_open + (length - 1) * scores.gap_extension

    if cigartuples:
        if cigartuples[0][0] != CIGAR_SOFTCLIP:
            score += scores.end_bonus
        if cigartuples[-1][0] != CIGAR_SOFTCLIP:
            score += scores.end_bonus

    return score


def filter_bam(alignment_file):
    for record in alignment_file:
        if not record.is_supplementary:
            yield record


def zip_longest_synthesize_unmapped(truth, predicted):
    p = None
    for t in truth:
        if p is None:
            p = next(predicted)
        tqn = t.query_name
        if t.is_paired and not t.is_read1:
            tqn += "/2"
        else:
            tqn += "/1"
        if p.query_name != tqn:
            synthetic = copy.copy(t)
            synthetic.is_unmapped = True
            synthetic.reference_start = 17
            synthetic.query_name = tqn

            yield t, synthetic
        else:
            yield t, p
            p = None


def get_iter_stats(truth, predicted, recompute_predicted_score=False, synthesize_unmapped=False):
    n = 0
    unaligned = 0
    nr_aligned = 0
    over_mapped = 0
    correct = 0
    correct_jaccard = 0.0
    correct_score = 0  # Same or better alignment score
    if synthesize_unmapped:
        iterator = zip_longest_synthesize_unmapped(filter_bam(truth), filter_bam(predicted))
    else:
        iterator = zip_longest(filter_bam(truth), filter_bam(predicted))

    for t, p in iterator:
        if t is None or p is None:
            raise ValueError("unequal number of records in the input files")

        if t.query_name != p.query_name:
            if p.query_name.endswith("/1") or p.query_name.endswith("/2"):
                p.query_name = p.query_name[:-2]
            if t.query_name != p.query_name:
                raise ValueError(f"query name mismatch: {t.query_name} != {p.query_name}")

        if t.is_paired != p.is_paired:
            PAIRED = {True: "paired", False: "single-end"}
            raise ValueError(
                f"Truth is {PAIRED[t.is_paired]} but predicted is {PAIRED[p.is_paired]}"
            )
        n += 1
        if t.is_unmapped:  # TODO and not p.is_unmapped:
            over_mapped += 1
            continue
        if not t.is_unmapped and p.is_unmapped:
            unaligned += 1
            continue
        if t.is_unmapped and p.is_unmapped:
            continue
        nr_aligned += 1

        is_correct = False
        if t.reference_name == p.reference_name:
            jacc = jaccard_overlap(
                p.reference_start, p.reference_end, t.reference_start, t.reference_end
            )
            correct_jaccard += jacc
            if overlap(
                p.reference_start, p.reference_end, t.reference_start, t.reference_end
            ):
                correct += 1
                is_correct = True

        if not is_correct:
            if recompute_predicted_score:
                predicted_score = recompute_alignment_score(p, Scores)
            else:
                predicted_score = p.get_tag("AS")
            truth_score = recompute_alignment_score(t, Scores)
            # print(f"true: {t.reference_name} {t.reference_start} {t.cigarstring} AS:{truth_score}  -- actual: {p.reference_name} {p.reference_start} {p.cigarstring} AS:{predicted_score}")
            correct_score += predicted_score >= truth_score

    aligned_percentage = 100 * nr_aligned / n
    accuracy = 100 * correct / n
    correct_score_percentage = 100 * (correct_score + correct) / n
    return (
        aligned_percentage,
        accuracy,
        over_mapped,
        100 * correct_jaccard / n,
        correct_score_percentage,
    )


def main(args):
    if args.predicted:
        with (
            AlignmentFile(args.truth) as truth,
            AlignmentFile(args.predicted) as predicted,
        ):
            if args.only_r1:
                truth = only_r1_iter(truth)
            if args.multiple_primary:
                if args.only_r1:
                    predicted = pick_random_primary_single_end_iter(predicted)
                else:
                    predicted = pick_random_primary_paired_end_iter(predicted)
            (
                percent_aligned,
                percent_correct,
                over_mapped,
                jacc,
                correct_score_percentage,
            ) = get_iter_stats(truth, predicted, args.recompute_score, synthesize_unmapped=args.synthesize_unmapped)
    elif args.predicted_paf:
        truth = read_alignments(args.truth, args.only_r1)
        predicted, mapped_to_multiple_pos = read_paf(args.predicted_paf)
        percent_aligned, percent_correct, over_mapped, jacc = get_stats(
            truth, predicted
        )
        correct_score = -1
        correct_score_percentage = -1

    print(
        percent_aligned,
        percent_correct,
        over_mapped,
        f"{jacc:.5}",
        correct_score_percentage,
        sep="\t",
    )
    # print("Percentage aligned: {}".format(round(percent_aligned, 3)))
    # print("Accuracy: {}".format(round(percent_correct, 3)))
    # print("Over-aligned (unmapped in grough truth file): {}".format(over_mapped))


if __name__ == "__main__":
    random.seed(0)
    parser = argparse.ArgumentParser(
        description="Calc identity",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--only-r1", default=False, action="store_true", help="Only use R1 reads from truth BAM")
    parser.add_argument("--recompute-score", default=False, action="store_true", help="Recompute score in *predicted* BAM. Default: Use score from AS tag")
    parser.add_argument("--multiple-primary", default=False, action="store_true", help="Allow multiple primary alignments (violates SAM specification) and pick one randomly")
    parser.add_argument("--synthesize-unmapped", default=False, action="store_true", help="If an alignment is missing from predicted, assume the read is unmapped")
    parser.add_argument("--truth", help="True SAM/BAM")
    parser.add_argument("--predicted", "--predicted_sam", help="Predicted SAM/BAM")
    parser.add_argument("--predicted_paf", help="Predicted PAF")
    parser.add_argument("--outfile", help="Path to file")

    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()

    main(args)
