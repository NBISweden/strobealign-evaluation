#!/usr/bin/env python3
"""
Compile a given commit of strobealign and store the binary in bin/
"""
from pathlib import Path
from tempfile import TemporaryDirectory
import shutil
import subprocess


def runc(*args, **kwargs):
    kwargs["check"] = True
    return subprocess.run(*args, **kwargs)



def compile_strobealign_if_missing(commit_hash, binary=None):
    """
    Compile the given strobealign commit and store the binary at
    *binary* (which must be a path). If *binary* is None, write to
    bin/strobealign-{short_commit_hash}. If the target binary already exists,
    do not compile.
    """
    if binary is None:
        short_commit_hash = make_short_hash(commit_hash)
        bindir = Path("bin")
        binary = bindir / f"strobealign-{short_commit_hash}"
    else:
        bindir = None

    if binary.exists():
        return binary, False

    bindir.mkdir(exist_ok=True)

    compile_strobealign(commit_hash, binary=binary)
    return binary, True


def compile_strobealign(commit_hash, binary):
    with TemporaryDirectory() as compiledir:
        runc(["git", "clone", "strobealign", compiledir])
        runc(["git", "checkout", "--detach", commit_hash], cwd=compiledir)

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


def make_short_hash(rev):
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
    parser.add_argument("hash", help="Commit hash")
    args = parser.parse_args()
    binary, was_built = compile_strobealign_if_missing(args.hash)


if __name__ == "__main__":
    main()
