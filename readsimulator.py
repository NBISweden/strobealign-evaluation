#!/usr/bin/env python3

# /// script
# dependencies = ["pyfaidx", "scipy"]
# ///

import sys
from argparse import ArgumentParser
import random
from pathlib import Path

from pyfaidx import Fasta
from scipy.stats import norm
import numpy


NEG = {
    "A": "CGT",
    "C": "AGT",
    "G": "ACT",
    "T": "ACG",
    "N": "N",
}

def mutate(s, error_rate):
    n = len(s)
    k = int(n * error_rate)
    positions = sorted(random.sample(range(len(s)), k=k))
    mutated = ""
    prev = 0
    md = ""
    for pos in positions:
        mutated += s[prev:pos]
        if pos - prev > 0:
            md += str(pos - prev)
        mutated_base = random.choice(NEG[s[pos]])
        mutated += mutated_base
        md += mutated_base
        prev = pos + 1
    mutated += s[prev:]
    if len(s) - prev > 0:
        md += str(len(s) - prev)

    assert len(s) == len(mutated)

    return mutated, md, k


def output_sam_records(name, seq1, seq2, contig, pos1, pos2, n_mismatches1, n_mismatches2, md1, md2):
    insert_size = pos2 - pos1 + len(seq2)
    print(
        name,
        99,  # PAIRED,PROPER_PAIR,MREVERSE,READ1
        contig,
        pos1 + 1,
        255,
        f"{len(seq1)}M",
        "=",
        pos2 + 1,
        insert_size,
        seq1,
        "H" * len(seq1),  # qual
        f"NM:i:{n_mismatches1}",
        f"MD:Z:{md1}",
        sep="\t",
    )
    print(
        name,
        147,  # PAIRED,PROPER_PAIR,REVERSE,READ2
        contig,
        pos2 + 1,
        255,
        f"{len(seq2)}M",
        "=",
        pos1 + 1,
        -insert_size,
        seq2,
        "H" * len(seq2),  # qual
        f"NM:i:{n_mismatches2}",
        f"MD:Z:{md2}",
        sep="\t",
    )


def output_sam_record(name, seq, contig, pos, n_mismatches, md):
    print(
        name,
        0,
        contig,
        pos + 1,
        255,
        f"{len(seq)}M",
        "*",
        0,
        0,
        seq,
        "H" * len(seq),  # qual
        f"NM:i:{n_mismatches}",
        f"MD:Z:{md}",
        sep="\t",
    )


def simulate_paired_end_reads(fasta, n, read_length, error_rate, mean_insert_size, stddev_insert_size):
    contig_names = list(fasta.keys())
    contig_lengths = [len(fr) for fr in fasta.values()]
    i = 0
    while i < n:
        contigs = random.choices(contig_names, weights=contig_lengths, k=10000)
        fragment_sizes = norm.rvs(loc=mean_insert_size, scale=stddev_insert_size, size=10000)
        for contig, fragment_size in zip(contigs, fragment_sizes):
            fragment_size = int(fragment_size)
            name = f"simulated.{i+1}"
            contig_length = len(fasta[contig])
            pos1 = random.randint(0, contig_length - fragment_size + 1)
            pos2 = pos1 + fragment_size - read_length
            assert pos2 >= 0
            seq1 = fasta[contig][pos1:pos1+read_length].seq.upper()
            seq2 = fasta[contig][pos2:pos2+read_length].seq.upper()
            if seq1.count("N") >= read_length / 10 or seq2.count("N") >= read_length / 10:
                continue
            if error_rate > 0:
                seq1, md1, n_mismatches1 = mutate(seq1, error_rate)
                seq2, md2, n_mismatches2 = mutate(seq2, error_rate)
            else:
                n_mismatches1 = n_mismatches2 = 0
                md1 = md2 = str(read_length)
            output_sam_records(name, seq1, seq2, contig, pos1, pos2, n_mismatches1, n_mismatches2, md1, md2)
            i += 1
            if i == n:
                break


def simulate_single_end_reads(fasta, n, read_length, error_rate):
    contig_names = list(fasta.keys())
    contig_lengths = [len(fr) for fr in fasta.values()]
    i = 0
    while i < n:
        contigs = random.choices(contig_names, weights=contig_lengths, k=10000)
        for contig in contigs:
            contig_length = len(fasta[contig])
            assert contig_length >= read_length
            name = f"simulated.{i+1}"
            pos = random.randint(0, len(fasta[contig]) - read_length + 1)
            seq = fasta[contig][pos:pos+read_length].seq.upper()
            if seq.count("N") >= read_length / 10:
                continue
            if error_rate > 0:
                seq, md, n_mismatches = mutate(seq, error_rate)
            else:
                n_mismatches = 0
                assert len(seq) == read_length
                md = str(read_length)
            output_sam_record(name, seq, contig, pos, n_mismatches, md)
            i += 1
            if i == n:
                break


def main():
    parser = ArgumentParser()
    parser.add_argument("-n", type=int, default=1_000_000, help="No. of read pairs to simulate")
    parser.add_argument("--error-rate", "-e", type=float, default=0, help="Error rate (only mismatches)")
    parser.add_argument("--se", dest="paired", action="store_false", default=True, help="Simulate single-end reads")
    parser.add_argument("--mean-insert-size", type=int, default=500)
    parser.add_argument("--stddev-insert-size", type=int, default=30)
    parser.add_argument("--read-length", type=int, default=300)
    parser.add_argument("--seed", type=int, default=0, help="Random seed")
    parser.add_argument("ref", type=Path, help="Reference FASTA")
    args = parser.parse_args()

    random.seed(args.seed)
    numpy.random.seed(args.seed)
    read_length = args.read_length

    with Fasta(args.ref) as fasta:
        records = {r.name: r for r in fasta if len(r) >= args.read_length}
        print(f"Filtered {len(fasta.keys()) - len(records)} contigs that are shorter than the read length", file=sys.stderr)
        if not records:
            raise ValueError("All contigs are shorter than the requested read length")
        print("@HD", "VN:1.4", sep="\t")
        for name, record in records.items():
            print("@SQ", f"SN:{name}", f"LN:{len(record)}", sep="\t")

        if args.paired:
            simulate_paired_end_reads(records, args.n, read_length, args.error_rate, args.mean_insert_size, args.stddev_insert_size)
        else:
            simulate_single_end_reads(records, args.n, read_length, args.error_rate)


if __name__ == "__main__":
    main()
