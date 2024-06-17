import os
import sys
import argparse
import pysam
from xopen import xopen


def read_sam(sam_file):
    SAM_file = pysam.AlignmentFile(sam_file, "r", check_sq=False)
    read_positions = {} # acc -> [ref_id, ref_start, refstop]

    for read in SAM_file.fetch(until_eof=True):
        if read.is_secondary:
            continue
        if read.flag == 0 or read.flag == 16: # single end
            if read.reference_end is None:
                read_positions[read.query_name] = False
            else:
                read_positions[read.query_name] = (read.reference_name, read.reference_start, read.reference_end)
        
        elif read.is_paired:
            q_name = read.query_name
            if read.is_read1 and not q_name.endswith("/1"):
                q_name += "/1"
        # elif read.flag == 99 or  read.flag == 83: # Paired end first
                # if not (read.flag == 99 or  read.flag == 83):
                #     print(read.query_name, read.flag)

        # elif read.flag == 147 or read.flag == 163: # Paired end second
            if read.is_read2 and not q_name.endswith("/2"):
                q_name += "/2"
                # if not (read.flag == 147 or read.flag == 163):
                #     print(read.query_name, read.flag)

            if read.is_unmapped:
                read_positions[q_name] = False
            else:
                read_positions[q_name] = (read.reference_name, read.reference_start, read.reference_end)
        
        elif read.is_unmapped:  # single and unmapped
            assert not read.is_paired
            read_positions[read.query_name] = False

    return read_positions        


def read_paf(paf_file):
    read_positions = {} # acc -> [ref_id, ref_start, refstop]
    mapped_to_multiple_pos = 0
    with xopen(paf_file) as paf:
        for line in paf:
            vals = line.split()
            read_acc, ref_name, reference_start, reference_end  = vals[0], vals[5], int(vals[7]), int(vals[8])

            if read_acc in read_positions:
                mapped_to_multiple_pos += 1
                continue
            else:
                read_positions[read_acc] = (ref_name, reference_start, reference_end)
    return read_positions, mapped_to_multiple_pos


def overlap(q_a, q_b, p_a, p_b):
    assert q_a <= q_b and p_a <= p_b
    # if (q_a == q_b) or (p_a == p_b):
    #     print("Cigar bug")
    return  (p_a <= q_a <= p_b) or (p_a <= q_b <= p_b) or (q_a <= p_a <= q_b) or (q_a <= p_b <= q_b)


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
    for read_acc in predicted:
        if not truth[read_acc]:
            over_mapped += 1
            continue
        if not predicted[read_acc]:
            unaligned += 1
            continue

        nr_aligned += 1

        pred_ref_id, pred_start, pred_stop = predicted[read_acc]
        true_ref_id, true_start, true_stop = truth[read_acc]

        # print(read_acc, pred_start, pred_stop, true_start, true_stop)
        if pred_ref_id == true_ref_id:
            if overlap(pred_start, pred_stop, true_start, true_stop):
                correct += 1
            correct_jaccard += jaccard_overlap(pred_start, pred_stop, true_start, true_stop)
    #         print(read_acc, pred_ref_id, pred_start, pred_stop, true_ref_id, true_start, true_stop )
    # print(correct)
    # print(nr_aligned)
    aligned_percentage = 100 * nr_aligned / nr_total
    accuracy = 100 * correct/nr_total
    return aligned_percentage, accuracy, over_mapped, 100 * correct_jaccard / nr_total


def main(args):
    truth = read_sam(args.truth)

    if args.predicted_sam:
        predicted = read_sam(args.predicted_sam)
    elif args.predicted_paf:
        predicted, mapped_to_multiple_pos = read_paf(args.predicted_paf)
        # print("Number of reads mapped to several positions (using first pos):", mapped_to_multiple_pos)

    percent_aligned, percent_correct, over_mapped, jacc = get_stats(truth, predicted)
    print(percent_aligned, percent_correct, over_mapped, f"{jacc:.5}", sep="\t")
    # print("Percentage aligned: {0}".format(round(percent_aligned, 3)))
    # print("Accuracy: {0}".format(round(percent_correct, 3)))
    # print("Over-aligned (unmapped in grough truth file): {0}".format(over_mapped))



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calc identity", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--truth', type=str, default=False, help='True coordinates (SAM)')
    parser.add_argument('--predicted_sam', type=str, default="", help='Perdicted coordinates (SAM/PAF)')
    parser.add_argument('--predicted_paf', type=str, default="", help='Perdicted coordinates (SAM/PAF)')
    parser.add_argument('--outfile', type=str, default=None, help='Path to file.')
    # parser.set_defaults(which='main')
    args = parser.parse_args()


    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()

    main(args)
