#!/usr/bin/env python3
"""
Compile a given commit of strobealign and store the binary in bin/ (or where
specified with -o)
"""
from pathlib import Path
from tempfile import TemporaryDirectory
import shutil
import subprocess

DEFAULT_REPOSITORY = "https://github.com/ksahlin/strobealign.git"


def runc(*args, **kwargs):
    kwargs["check"] = True
    return subprocess.run(*args, **kwargs)



def compile_strobealign_if_missing(commit_hash, binary=None, repository=None):
    """
    Compile the given strobealign commit and store the binary at
    *binary* (which must be a path). If *binary* is None, write to
    bin/strobealign-{short_commit_hash}. If the target binary already exists,
    do not compile.
    """
    if repository is None:
        repository = DEFAULT_REPOSITORY
    if binary is None:
        #short_commit_hash = make_short_hash(commit_hash)
        bindir = Path("bin")
        bindir.mkdir(exist_ok=True)
        binary = bindir / f"strobealign-{commit_hash}"

    if binary.exists():
        return binary, False

    success = False
    compile_strobealign(commit_hash, binary=binary, repository=repository)

    return binary, True


class CouldNotCheckout(Exception):
    pass


def compile_strobealign(commit_hash, binary, repository=DEFAULT_REPOSITORY):
    with TemporaryDirectory() as compiledir:
        runc(["git", "clone", repository, compiledir])
        if subprocess.run(["git", "checkout", "--detach", commit_hash], cwd=compiledir).returncode != 0:
            raise CouldNotCheckout()

        runc(
            [
                "cmake",
                compiledir,
                "-DCMAKE_RULE_MESSAGES=OFF",
                "--log-level=NOTICE",
                '-DCMAKE_C_FLAGS="-march=native"',
                '-DCMAKE_CXX_FLAGS="-march=native"',
                "-B",
                Path(compiledir) / "build",
            ]
        )
        runc(
            [
                "cmake",
                "--build",
                f"{compiledir}/build",
                "-j",
                "8",
                "--target",
                "strobealign",
            ]
        )
        shutil.move(Path(compiledir) / "build" / "strobealign", binary)


def _currently_unused_make_short_hash(rev):
    result = runc(
        ["git", "rev-parse", "--short", f"{rev}^0"],
        cwd="strobealign",
        capture_output=True,
        text=True,
    )
    short_hash = result.stdout.strip()
    return short_hash


def main():
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("-o", "--output", type=Path, help="Output file name")
    parser.add_argument("-r", "--repository", default=DEFAULT_REPOSITORY, help="Git repository from which to clone. Default: %(default)s")
    parser.add_argument("hash", help="Commit hash")

    args = parser.parse_args()
    binary, was_built = compile_strobealign_if_missing(args.hash, args.output, repository=args.repository)


if __name__ == "__main__":
    main()
