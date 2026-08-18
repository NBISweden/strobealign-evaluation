#!/usr/bin/env python3
"""
Benchmark two strobealign commits against each other
"""

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

GENOMES = ["fruitfly", "maize" "CHM13", "chrY", "rye"]
DATASETS = ["sim4", "sim6", "sim0"]
READ_LENGTHS = [50, 100, 150, 200, 300, 500, 1000, 5000, 10000]


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
        "commit1", help="First commit (can also be a tag or branch name)"
    )
    parser.add_argument(
        "commit2",
        default=None,
        nargs="?",
        help="Second commit. If not provided, compare merge base of main and commit1 to commit1",
    )
    args = parser.parse_args()
    logging.basicConfig(format="%(message)s", level="DEBUG" if args.debug else "INFO")
    del args.debug
    run(**vars(args))


def run(commit1: str, commit2: str | None, config: Path):
    config = read_config(config)
    logger.debug("Configuration: %s", config)

    if commit2 is None:
        commit1, commit2 = get_merge_base(Path(config.repository), commit1), commit1

    logger.info("Comparing commit %s to %s", commit1, commit2)

    # Compile
    programs = [
        make_program(commit, binary_dir=config.binaries, repository=config.repository)
        for commit in [commit1, commit2]
    ]

    modes = [Mode.MAP, Mode.ALIGN]
    libraries = make_libraries()

    print(
        "dataset",
        f"{'library':14}",
        "ends",
        "mode",
        "accuracy",
        "accuracy_diff",
        "speedup",
        "avg_accuracy_diff",
        "avg_speedup",
        sep="\t",
    )
    accuracy_diff_sum = 0.0
    speedup_sum = 0.0
    n = 0
    for mode, library in itertools.product(modes, libraries):
        accuracies, runtimes = [], []
        for program in programs:
            accuracy, runtime = run_one_library(program, config, library, mode)
            accuracies.append(accuracy)
            runtimes.append(runtime)

        n += 1
        accuracy_diff_sum += accuracies[1].correct - accuracies[0].correct
        speedup = runtimes[0].mapping_time / runtimes[1].mapping_time
        speedup_sum += speedup
        library_name = f"{library.genome}-{library.read_length}"
        print(
            library.dataset,
            f"{library_name:14}",
            library.ends,
            "map  " if mode == Mode.MAP else "align",
            f"{accuracies[0].correct:9.3f}",
            f"{accuracies[1].correct - accuracies[0].correct:9.4f}",
            f"{speedup:6.2f}",
            f"{accuracy_diff_sum / n:+.6f}",
            f"{speedup_sum / n:6.3f}",
            sep="\t",
        )

    print("Average accuracy change:", f"{accuracy_diff_sum / n:+.6f}", file=sys.stderr)

    # TODO remove output file


def make_libraries() -> list[Library]:
    libraries = []
    for dataset in DATASETS:
        for genome in GENOMES:
            for read_length in READ_LENGTHS:
                # TODO
                if not (
                    Path("datasets") / dataset / f"{genome}-{read_length}" / "truth.bam"
                ).exists():
                    continue
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
