#!/usr/bin/env python3
"""
Benchmark two strobealign commits against each other

Alternatively, benchmark different command-line arguments against each other
"""

from enum import Enum
import itertools
import subprocess

# TODO
# Output is a table with the following columns:
# - dataset:
# - library
# - ends (se for single-end, pe for paired-end)
# - mode (map or align)
# - accuracy_diff: difference of accuracies, in percentage points
# - runtime_diff: difference of runtimes, excluding indexing, in seconds
# """

import sys
from argparse import ArgumentParser
from pathlib import Path
import logging

import pandas as pd

from .run1 import (
    Library,
    Mode,
    Ends,
    read_config,
    make_program,
    Job,
)
from .plots import plot as saeval_plot, MEASUREMENT_TYPES

logger = logging.getLogger(__name__)

GENOMES = ["blueflax", "fruitfly", "maize", "CHM13", "rye", "chrY"]
DATASETS = ["sim6", "sim0", "sim4", "sim6"]
READ_LENGTHS = [50, 100, 150, 200, 300, 500, 1000, 5000, 10000]
LONG_READ_LENGTH = 1000  # no paired-end


class AccuracyMethod(Enum):
    overlap = "overlap"
    jaccard = "jaccard"
    score = "score"

    def __str__(self):
        return self.value


def main():
    parser = ArgumentParser()
    parser.add_argument("--debug", action="store_true", default=False)
    parser.add_argument(
        "--config",
        type=Path,
        help="Path to saeval configuration file. "
        "Default: $XDG_CONFIG_HOME/saeval/config.toml or ~/.config/saeval/config.toml",
    )
    parser.add_argument(
        "--markdown",
        default=False,
        action="store_true",
        help="Output table in Markdown format",
    )
    parser.add_argument("--plot", type=Path, metavar="PDF", help="Plot results to PDF")
    parser.add_argument(
        "--accuracy",
        type=AccuracyMethod,
        default=AccuracyMethod.overlap,
        choices=list(AccuracyMethod),
        help="Accuracy computation method. Default: overlap",
    )
    parser.add_argument("--mode", dest="modes", choices=["map", "align"])
    parser.add_argument("--ends", dest="ends", choices=["se", "pe"])
    parser.add_argument(
        "--args1",
        help="Arguments to pass to first strobealign commit. "
        'Write as --args1= (with an "=") like this: --args1="--arg1=value --arg2=value"',
    )
    parser.add_argument(
        "--args2", help="Arguments to pass to second strobealign commit"
    )
    parser.add_argument(
        "commit1", help="First commit (can also be a tag or branch name)"
    )
    parser.add_argument(
        "commit2",
        default=None,
        nargs="?",
        help="Second commit. If not provided, compare merge base of main and commit1 to commit1",
    )
    args = parser.parse_args()
    if args.modes == "map":
        modes = [Mode.MAP]
    elif args.modes == "align":
        modes = [Mode.ALIGN]
    else:
        modes = [Mode.MAP, Mode.ALIGN]
    args.modes = modes

    if args.ends == "se":
        ends = [Ends.SE]
    elif args.ends == "pe":
        ends = [Ends.PE]
    else:
        ends = [Ends.SE, Ends.PE]
    args.ends = ends

    logging.basicConfig(format="%(message)s", level="DEBUG" if args.debug else "INFO")
    del args.debug
    run(**vars(args))


def run(
    commit1: str,
    commit2: str | None,
    config: Path,
    args1: str | None,
    args2: str | None,
    accuracy: AccuracyMethod,
    modes: list[Mode],
    ends: list[Ends],
    markdown: bool,
    plot: Path | None,
):
    config = read_config(config)
    logger.debug("Configuration: %s", config)

    if commit2 is None:
        commit1, commit2 = get_merge_base(Path(config.repository), commit1), commit1

    # Compile
    programs = [
        make_program(commit, binary_dir=config.binaries, repository=config.repository)
        for commit in [commit1, commit2]
    ]
    logger.info("Comparing commit %s to %s\n", commit1, commit2)
    if args1 or args2:
        logger.info(
            "Command-line arguments: '%s' vs '%s'",
            args1 if args1 else "",
            args2 if args2 else "",
        )
        args1 = args1.split() if args1 else []
        args2 = args2.split() if args2 else []

    libraries = make_libraries(ends)

    if accuracy is AccuracyMethod.overlap:
        col = "accuracy"
        acc = lambda a: a.correct
    elif accuracy is AccuracyMethod.score:
        col = "saccuracy"
        acc = lambda a: a.score_correct
    elif accuracy is AccuracyMethod.jaccard:
        col = "jaccuracy"
        acc = lambda a: a.jaccard_correct

    if markdown:
        sep = " | "
    else:
        sep = "\t"

    cols = [
        "mode ",
        f"{'library':22}",
        f"{col} (%)",
        f"∆{col}",
        "speedup",
    ]
    if not markdown:
        cols += [
            f"avg_∆{col}",
            "avg_speedup",
        ]
    print(*cols, sep=sep)
    if markdown:
        print("-|-|-:|-:|-:")
    accuracy_diff_sum = 0.0
    speedup_sum = 0.0
    n = 0
    rows = []
    for mode, library in itertools.product(modes, libraries):
        library_name = f"{library.genome}-{library.read_length}"
        print(
            "map  " if mode == Mode.MAP else "align",
            f"{library.ends}-{library.dataset}-{library_name}".ljust(22),
            sep=sep,
            end=sep,
        )
        sys.stdout.flush()
        job1 = Job(programs[0], args1, config, library, mode)
        accuracy1, runtime1 = job1.run()
        print(f"{acc(accuracy1):9.3f}", end=sep)
        sys.stdout.flush()

        job2 = Job(programs[1], args2, config, library, mode)
        accuracy2, runtime2 = job2.run()
        n += 1
        accuracy_diff_sum += acc(accuracy2) - acc(accuracy1)
        speedup = runtime1.mapping_time / runtime2.mapping_time
        speedup_sum += speedup
        cols = [
            f"{acc(accuracy2) - acc(accuracy1):9.4f}",
            f"{speedup:6.2f}×",
        ]
        if not markdown:
            cols += [
                f"{accuracy_diff_sum / n:+.6f}",
                f"{speedup_sum / n:6.3f}",
            ]
        print(*cols, sep=sep)
        rows.append(to_dict(job1, accuracy1, runtime1))
        rows.append(to_dict(job2, accuracy2, runtime2))

    print()
    print(f"Average accuracy change: {accuracy_diff_sum / n:+.6f}", file=sys.stderr)
    print(f"Average speedup: {speedup_sum / n:6.3f}", file=sys.stderr)

    if plot:
        import pandas as pd

        df = pd.DataFrame(rows)
        print(df.to_csv(index=False))
        plot_ends_to_pdf(df, plot)

    # TODO remove output file


def to_dict(job, accuracy, runtime):
    library = job.library
    d = dict()
    d["tool"] = job.tool_tag()
    d["mode"] = str(job.mode)
    d["dataset"] = library.dataset
    d["genome"] = library.genome
    d["ends"] = str(library.ends)
    d["read_length"] = library.read_length
    d["read_count"] = runtime.reads
    d["aligned"] = accuracy.aligned
    d["overaligned"] = accuracy.overmapped
    d["accuracy"] = accuracy.correct
    d["jaccuracy"] = accuracy.jaccard_correct
    d["saccuracy"] = accuracy.score_correct
    d["time"] = runtime.mapping_time
    d["memory"] = round(runtime.memory_mib / 1024, 3)

    return d


def make_libraries(ends_list: list[Ends]) -> list[Library]:
    libraries = []
    all = itertools.product(ends_list, DATASETS, GENOMES, READ_LENGTHS)
    for ends, dataset, genome, read_length in all:
        if ends == Ends.PE and read_length >= LONG_READ_LENGTH:
            continue
        libraries.append(
            Library(
                dataset=dataset,
                genome=genome,
                read_length=read_length,
                ends=ends,
            )
        )

    return libraries


def get_merge_base(repo: Path, commit: str) -> str:
    result = subprocess.run(
        ["git", "merge-base", commit, "main"],
        cwd=repo,
        capture_output=True,
        text=True,
        check=True,
    )
    base = result.stdout.strip()

    return base


def plot_ends_to_pdf(df, path):
    import seaborn as sns
    import matplotlib
    from matplotlib.backends.backend_pdf import PdfPages

    matplotlib.use("PDF")

    sns.set(font_scale=1.4)
    sns.set_style("whitegrid")

    df["time"] = df["time"] * 1e6 / df["read_count"]

    modes = sorted(set(df["mode"]))
    datasets = sorted(set(df["dataset"]))
    read_lengths = sorted(set(df["read_length"]))
    tools = {tag: tag.removeprefix("strobealign-") for tag in sorted(set(df["tool"]))}

    palette = {
        df["tool"][0]: "tab:green",
        df["tool"][1]: "black",
    }

    with PdfPages(path) as pdf:
        for ends, table in df.groupby("ends"):
            title = "Single-end reads" if ends == "se" else "Paired-end reads"
            for y, label, logscale, in_pdf in MEASUREMENT_TYPES:
                fig = saeval_plot(
                    table,
                    palette,
                    tools,
                    modes,
                    datasets,
                    read_lengths,
                    y=y,
                    logscale=logscale,
                    xlogscale=True,
                    row="genome",
                    label=label,
                    title=title,
                    linewidth=2,
                    solid=False,
                    legend="right",
                )
                if in_pdf:
                    pdf.savefig(bbox_inches="tight")


if __name__ == "__main__":
    main()
