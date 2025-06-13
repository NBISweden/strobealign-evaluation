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


def output_sam_records(name, seq1, seq2, contig, pos1, pos2):
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
        "NM:i:0",
        f"MD:Z:{len(seq1)}",
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
        "NM:i:0",
        f"MD:Z:{len(seq2)}",
        sep="\t",
    )


def output_sam_record(name, seq, contig, pos):
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
        "NM:i:0",
        f"MD:Z:{len(seq)}",
        sep="\t",
    )


def simulate_paired_end_reads(fasta, n, read_length, mean_insert_size, stddev_insert_size):
    contig_names = list(fasta.keys())
    contig_lengths = [len(fr) for fr in fasta]
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
            output_sam_records(name, seq1, seq2, contig, pos1, pos2)
            i += 1
            if i == n:
                break


def simulate_single_end_reads(fasta, n, read_length):
    contig_names = list(fasta.keys())
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
            output_sam_record(name, seq, contig, pos)
            i += 1
            if i == n:
                break


def main():
    parser = ArgumentParser()
    parser.add_argument("-n", type=int, default=1_000_000, help="No. of read pairs to simulate")
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
            simulate_paired_end_reads(records, args.n, read_length, args.mean_insert_size, args.stddev_insert_size)
        else:
            simulate_single_end_reads(records, args.n, read_length)


if __name__ == "__main__":
    main()
