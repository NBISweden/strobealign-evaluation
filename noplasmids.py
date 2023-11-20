#!/usr/bin/env python3
import sys
import dnaio

def print_fasta(record, line_length=80):
    print(f">{record.name}")
    for i in range(0, len(record.sequence), line_length):
        print(record.sequence[i:i+line_length])


for path in sys.argv[1:]:
    with dnaio.open(path) as inf, dnaio.open(sys.stdout.buffer, mode="w", fileformat="fasta") as outf:
        for record in inf:
            if "chromosome" in record.name:
                print_fasta(record)
