#!/usr/bin/env python3
"""
Combine two BAM files from two strobealign runs with the aim of
improving mapping rate and accuracy
"""
from argparse import ArgumentParser
from contextlib import ExitStack

from pysam import AlignmentFile


def main():
    parser = ArgumentParser()
    parser.add_argument("-s", "--single-end", action="store_true")
    parser.add_argument("--output")
    parser.add_argument("bam1")
    parser.add_argument("bam2")
    args = parser.parse_args()

    with ExitStack() as stack:
        bam1 = stack.enter_context(AlignmentFile(args.bam1))
        bam2 = stack.enter_context(AlignmentFile(args.bam2))
        out = stack.enter_context(AlignmentFile(args.output, "wb", template=bam1))

        if args.single_end:
            combine_single_end(bam1, bam2, out)
        else:
            combine_paired_end(bam1, bam2, out)


def combine_single_end(bam1, bam2, out):
    for left, right in zip(bam1, bam2):
        which = 0
        if left.is_unmapped:
            which = 1
        elif right.is_unmapped:
            which = 0
        elif left.get_tag("AS") < right.get_tag("AS"):
            which = 1

        out.write((left, right)[which])


def combine_paired_end(bam1, bam2, out):
    by_proper = 0
    by_unmapped = 0
    by_score = 0
    # by_contig = 0

    for (left1, left2), (right1, right2) in zip(pairs(bam1), pairs(bam2)):
        left_mapped = 2 - left1.is_unmapped - left2.is_unmapped
        right_mapped = 2 - right1.is_unmapped - right2.is_unmapped

        which = 0
        try:
            left_score = left1.get_tag("AS") + left2.get_tag("AS")
            right_score = right1.get_tag("AS") + right2.get_tag("AS")
            if right_score > left_score:
                by_score += 1
                which = 1
        except KeyError:
            pass
        if not left1.is_proper_pair and right1.is_proper_pair:
            by_proper += 1
            which = 1
        elif right_mapped > left_mapped:
            by_unmapped += 1
            which = 1
        #elif left1.reference_name != left2.reference_name and right1.reference_name == right2.reference_name:
        #    by_contig += 1
        #    which = 1

        out.write((left1, right1)[which])
        out.write((left2, right2)[which])
    # print("Reasons:")
    # print("- proper pair:", by_proper)
    # print("- unmapped:", by_unmapped)
    # print("- score:", by_score)
    print(by_proper, by_unmapped, by_score, sep="\t")


def pairs(bam):
    it = iter(bam)
    for r1, r2 in zip(it, it):
        assert r1.is_read1
        assert r2.is_read2
        yield r1, r2


if __name__ == "__main__":
    main()
