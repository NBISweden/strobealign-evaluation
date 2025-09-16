from dataclasses import dataclass


@dataclass
class Runtimes:
    user_time: float
    memory_mib: float
    indexing_time: float
    mapping_time: float  # This is what we usually report
    threads: int


def parse_log(log) -> Runtimes:
    # user_time = memory_mib = indexing_time = mapping_time = threads = None
    with open(log) as f:
        for line in f:
            if "User time (seconds)" in line:
                user_time = float(line.split()[-1])
            elif "Maximum resident" in line:
                memory_mib = float(line.split()[-1]) / 1024
            elif "Total time indexing:" in line:
                indexing_time = float(line.split()[3])
            elif "Total time mapping:" in line:
                mapping_time = float(line.split()[3])
            elif line.startswith("Threads:"):
                threads = int(line.split()[1])
            elif line.startswith("Mapping threads:"):
                threads = int(line.split()[2])
    return Runtimes(
        user_time=user_time,
        memory_mib=memory_mib,
        indexing_time=indexing_time,
        mapping_time=mapping_time,
        threads=threads,
    )

