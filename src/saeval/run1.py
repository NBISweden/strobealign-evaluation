#!/usr/bin/env python3
"""
Run strobealign on a single library and output accuracy and runtime as CSV.
"""

from argparse import ArgumentParser
from contextlib import ExitStack
from dataclasses import dataclass
from enum import Enum, StrEnum
import logging
import os
from pathlib import Path
import subprocess
import tomllib

from .compile import DEFAULT_REPOSITORY, compile_strobealign_if_missing, make_short_hash
from .get_accuracy import measure_accuracy, AccuracyPercentages
from .logs import parse_log

logger = logging.getLogger(__name__)


DEFAULT_THREADS = 8


class Mode(Enum):
    ALIGN = 0
    MAP = 1


class Ends(StrEnum):
    SE = "se"
    PE = "pe"


@dataclass
class Library:
    ends: Ends
    dataset: str
    genome: str
    read_length: int


@dataclass
class Program:
    commit: str
    binary: Path


@dataclass
class Configuration:
    repository: str
    genomes: Path
    datasets: Path
    binaries: Path
    results: Path
    threads: int


def main():
    parser = ArgumentParser()
    parser.add_argument("--debug", action="store_true", default=False)
    parser.add_argument(
        "--config",
        type=Path,
        help="Path to saeval configuration file. Default: $XDG_CONFIG_HOME/saeval/config.toml or ~/.config/saeval/config.toml",
    )
    # parser.add_argument("--ends", choices=("se", "pe"), default="se")
    parser.add_argument("--mode", choices=("align", "map"), default="align")
    parser.add_argument("--commit", help="strobealign commit", default="main")
    parser.add_argument(
        "library",
        help="Library to run in the form <ends>-<dataset>-<genome>-<readlength> (e.g. se-sim6-fruitfly-100)",
    )
    args = parser.parse_args()
    logging.basicConfig(format="%(message)s", level="DEBUG" if args.debug else "INFO")

    config = read_config(args.config)
    logger.debug("Configuration: %s", config)

    # Compile
    if args.commit is None:
        parser.error("Commit required (use -c)")
    program = make_program(
        args.commit, binary_dir=config.binaries, repository=config.repository
    )

    library = parse_library_tag(args.library)
    mode = Mode.ALIGN if args.mode == "align" else Mode.MAP

    logger.debug("Library: %s", library)
    read_count = read_read_count(config, library)

    args = []
    accuracy, runtime = Job(program, args, config, library, mode).run()

    logger.debug("Accuracy: %s", accuracy)
    logger.debug("Runtime: %s", runtime)

    # Example:
    # strobealign-16bf9274,align,sim0,fruitfly,se,50,1000000,99.9863,0,89.0574,89.0537,99.9396,1.58,0.821148
    print(
        f"strobealign-{program.commit}",
        "align" if mode == Mode.ALIGN else "map",
        library.dataset,
        library.genome,
        library.ends,
        library.read_length,
        read_count,
        accuracy.aligned,
        accuracy.overmapped,
        accuracy.correct,
        accuracy.jaccard_correct,
        accuracy.score_correct,
        runtime.mapping_time,
        f"{runtime.memory_mib / 1024:.6}",
        sep=",",
    )


def parse_library_tag(tag: str) -> Library:
    fields = tag.split("-")
    if len(fields) != 4:
        raise ValueError("Library tag needs to have four dash-separated fields")
    ends, dataset, genome, read_length = fields
    ends = Ends(ends)
    read_length = int(read_length)

    return Library(ends=ends, dataset=dataset, genome=genome, read_length=read_length)


def read_config(path: Path | None = None) -> Configuration:
    """
    Read saeval configuration.

    If path is not provided, read from
    $XDG_CONFIG_HOME/.config/saeval/config.toml
    """
    if path is None:
        config_home = os.getenv("XDG_CONFIG_HOME", os.path.expanduser("~/.config"))
        path = Path(config_home) / "saeval" / "config.toml"
    logger.debug("Reading configuration from %s", path)
    with open(path, "rb") as f:
        parsed = tomllib.load(f)

    repo = parsed.get("repository", DEFAULT_REPOSITORY)
    if "~" in repo and "://" not in repo:
        repo = str(Path(repo).expanduser())

    return Configuration(
        repository=repo,
        genomes=Path(parsed["genomes"]).expanduser(),
        datasets=Path(parsed["datasets"]).expanduser(),
        binaries=Path(parsed["binaries"]).expanduser(),
        results=Path(parsed["results"]).expanduser(),
        threads=parsed.get("threads", 8),
    )


def make_program(commit, binary_dir, repository) -> Program:
    if "://" not in repository:
        # `git rev-parse` does not work with a remote repository.
        short_hash = make_short_hash(Path(repository), commit)
    else:
        short_hash = commit
    (binary_dir / "strobealign").mkdir(exist_ok=True, parents=True)
    binary = binary_dir / "strobealign" / short_hash

    _, was_built = compile_strobealign_if_missing(
        short_hash, binary=binary, repository=repository
    )
    logger.debug("%s", f"Built {commit} ({short_hash})")
    return Program(binary=binary, commit=short_hash)


def read_read_count(config: Configuration, library: Library) -> int:
    libdir = (
        config.datasets / library.dataset / f"{library.genome}-{library.read_length}"
    )
    read_count = int((libdir / f"{library.ends}.readcount.txt").read_text())

    return read_count


class Job:
    """
    A single run of strobealign
    """
    def __init__(self,
        program: Program,
        args: list[str] | None,
        config: Configuration,
        library: Library,
        mode: Mode,
    ):
        self.program = program
        self.args = args
        self.config = config
        self.library = library
        self.mode = mode

    def run(self) -> (float, float):
        """
        Return the results for this job (accuracy and runtime), running the job
        if necessary.
        """
        outdir = self._outdir()

        logger.debug("Output directory: %s", outdir)
        ends = self.library.ends
        if self.mode == Mode.MAP:
            log_txt = outdir / f"{ends}.paf.log"
            outpath = outdir / f"{ends}.paf.gz"
            accuracy_txt = outdir / f"{ends}.paf.accuracy.txt"
        else:
            log_txt = outdir / f"{ends}.bam.log"
            outpath = outdir / f"{ends}.bam"
            accuracy_txt = outdir / f"{ends}.bam.accuracy.txt"

        if not log_txt.exists() or not accuracy_txt.exists():
            libdir = (
                self.config.datasets / self.library.dataset / f"{self.library.genome}-{self.library.read_length}"
            )
            ref = self.config.genomes / f"{self.library.genome}.fa"
            reads = [libdir / "1.fastq.gz"]
            if self.library.ends == Ends.PE:
                reads.append(libdir / "2.fastq.gz")
            truth = libdir / "truth.bam"
            if not truth.exists():
                raise FileNotFoundError(f"Truth {truth} not available")
            if not ref.exists():
                raise FileNotFoundError(f"Reference FASTA {ref} not found")

            outdir.mkdir(parents=True, exist_ok=True)
            if not outpath.exists() or not log_txt.exists():
                run_strobealign(
                    self.program.binary,
                    ref,
                    reads,
                    outpath,
                    log_txt,
                    self.mode,
                    self.config.threads,
                    extra_args=self.args,
                )
            assert outpath.exists() and log_txt.exists()
            logger.debug("Getting accuracy")
            accuracy = measure_accuracy(
                truth=truth, predicted=outpath, skip_r2=self.library.ends == Ends.SE
            )
            with open(accuracy_txt, mode="w") as f:
                print(*accuracy.percentages().row(), sep="\t", file=f)

        logger.debug("Reading %s", accuracy_txt)
        accuracy = parse_accuracy(accuracy_txt)
        runtime = parse_log(log_txt)

        return accuracy, runtime


    def _outdir(self) -> Path:
        libname = f"{self.library.dataset}-{self.library.genome}-{self.library.read_length}"
        key = self.program.commit
        if self.args:
            # Normalize
            key += "-" + "".join(arg.replace("-", "").replace("=", "") for arg in self.args)

        return self.config.results / f"strobealign-{key}" / libname


def run_strobealign(
    binary: Path,
    ref: Path,
    reads: list[Path],
    outpath: Path,
    log: Path,
    mode: Mode,
    threads: int = 8,
    extra_args: list[str] = None,
):
    tmp_log = log.with_suffix(".log.tmp")
    args = ["/usr/bin/time", "-v", binary, "-v"]
    if mode == Mode.MAP:
        args += ["-x"]
    if extra_args:
        args += extra_args
    args += [
        "-t",
        str(threads),
        str(ref),
        *(str(r) for r in reads),
    ]

    # Fill disk caches
    for path in binary, ref, *reads:
        subprocess.run(["cat", path], stdout=subprocess.DEVNULL, check=True)

    tmp_outpath = outpath.with_name(outpath.stem + ".tmp" + outpath.suffix)
    with ExitStack() as stack:
        logf = stack.enter_context(open(tmp_log, "wb"))
        logger.debug("Log to: %s", tmp_log)
        logger.debug("Running: %s", " ".join(str(arg) for arg in args))
        strobealign_process = subprocess.Popen(
            args, stdout=subprocess.PIPE, stderr=logf
        )

        if mode == Mode.MAP:
            outfile = stack.enter_context(open(tmp_outpath, "wb"))
            compress_process = subprocess.Popen(
                ["igzip"], stdin=strobealign_process.stdout, stdout=outfile
            )
        else:
            compress_process = subprocess.Popen(
                ["samtools", "view", "--no-PG", "-o", tmp_outpath],
                stdin=strobealign_process.stdout,
            )

        strobealign_process.wait()
        compress_process.wait()
        if strobealign_process.returncode != 0:
            raise ValueError()

    tmp_outpath.rename(outpath)
    tmp_log.rename(log)


def parse_accuracy(path) -> AccuracyPercentages:
    text = path.read_text()
    fields = text[:-1].split("\t")
    return AccuracyPercentages(
        aligned=float(fields[0]),
        correct=float(fields[1]),
        overmapped=int(fields[2]),
        jaccard_correct=float(fields[3]),
        score_correct=None if fields[4] == "" else float(fields[4]),
    )


if __name__ == "__main__":
    main()
