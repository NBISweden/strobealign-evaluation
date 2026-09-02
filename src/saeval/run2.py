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

from .run1 import (
    Library,
    Mode,
    Ends,
    read_config,
    make_program,
    run_one_library,
)

logger = logging.getLogger(__name__)

GENOMES = ["fruitfly", "maize", "CHM13", "chrY", "rye"]
DATASETS = ["sim4", "sim6", "sim0"]
READ_LENGTHS = [50, 100, 150, 200, 300, 500, 1000, 5000, 10000]


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
    parser.add_argument(
        "--accuracy",
        type=AccuracyMethod,
        default=AccuracyMethod.overlap,
        choices=list(AccuracyMethod),
        help="Accuracy computation method. Default: overlap",
    )
    parser.add_argument("--mode", dest="modes", choices=["map", "align"])
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
    markdown: bool,
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

    libraries = make_libraries()

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
    for mode, library in itertools.product(modes, libraries):
        accuracies, runtimes = [], []
        for program, args in zip(programs, [args1, args2]):
            accuracy, runtime = run_one_library(program, args, config, library, mode)
            accuracies.append(accuracy)
            runtimes.append(runtime)

        n += 1
        accuracy_diff_sum += acc(accuracies[1]) - acc(accuracies[0])
        speedup = runtimes[0].mapping_time / runtimes[1].mapping_time
        speedup_sum += speedup
        library_name = f"{library.genome}-{library.read_length}"

        cols = [
            "map  " if mode == Mode.MAP else "align",
            f"{library.ends}-{library.dataset}-{library_name}".ljust(22),
            f"{acc(accuracies[0]):9.3f}",
            f"{acc(accuracies[1]) - acc(accuracies[0]):9.4f}",
            f"{speedup:6.2f}×",
        ]
        if not markdown:
            cols += [
                f"{accuracy_diff_sum / n:+.6f}",
                f"{speedup_sum / n:6.3f}",
            ]
        print(*cols, sep=sep)

    print()
    print("Average accuracy change:", f"{accuracy_diff_sum / n:+.6f}", file=sys.stderr)

    # TODO remove output file


def make_libraries() -> list[Library]:
    libraries = []
    all = itertools.product(DATASETS, GENOMES, READ_LENGTHS)
    for dataset, genome, read_length in all:
        libraries.append(
            Library(
                dataset=dataset,
                genome=genome,
                read_length=read_length,
                ends=Ends.SE,
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


if __name__ == "__main__":
    main()
