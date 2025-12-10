#!/usr/bin/env python3
"""
Split chromosomes larger than 1 Gbp into multiple records
"""
import sys

import dnaio

max_chunk_size = 1_000_000_000

with (
    dnaio.open(sys.argv[1]) as infile,
    dnaio.FastaWriter(sys.stdout.buffer, line_length=80) as outfile
):
    for record in infile:
        if len(record) <= max_chunk_size:
            outfile.write(record)
            continue

        chunk_starts = range(0, len(record), max_chunk_size)
        n_chunks = len(chunk_starts)
        for i, start in enumerate(chunk_starts, 1):
            chunk = record[start:start+max_chunk_size]
            chunk.name = f"{record.name}_{i}of{n_chunks} chunk len={len(chunk)}"
            outfile.write(chunk)
