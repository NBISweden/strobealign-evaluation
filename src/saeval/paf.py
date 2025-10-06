from dataclasses import dataclass
from .get_accuracy import ReferenceInterval


@dataclass
class PafRecord:
    query_name: str
    query_length: int
    query_start: int
    query_end: int
    strand: str # ??
    ref_interval: ReferenceInterval
    n_matches: int
    block_length: int
    mapq: int
    tags: dict


def parse_paf(file):
    for line in file:
        fields = line.split("\t")
        query_name, query_length, query_start, query_end, strand, ref_name, ref_length, ref_start, ref_end, n_matches, block_length, mapq, *tags = fields
        yield PafRecord(
            query_name=query_name,
            query_length=int(query_length),
            query_start=int(query_start),
            query_end=int(query_end),
            strand=strand,
            ref_interval=ReferenceInterval(ref_name, int(ref_start), int(ref_end)),
            n_matches=int(n_matches),
            block_length=int(block_length),
            mapq=int(mapq),
            tags=tags,
        )
