#!/usr/bin/env python3
"""
Compile a given commit of strobealign and store the binary in bin/
"""
from pathlib import Path
from tempfile import TemporaryDirectory
import subprocess


def runc(*args, **kwargs):
    kwargs["check"] = True
    return subprocess.run(*args, **kwargs)


def compile_strobealign(commit_hash):
    short_commit_hash = make_short_hash(commit_hash)

    bindir = Path("bin")
    binary = bindir / f"strobealign-{short_commit_hash}"

    if binary.exists():
        return binary, False

    Path("builds").mkdir(exist_ok=True)
    with TemporaryDirectory(dir="builds") as compiledir:
        runc(["git", "clone", "strobealign", compiledir])
        runc(["git", "checkout", "--detach", short_commit_hash], cwd=compiledir)

        runc(
            [
                "cmake",
                compiledir,
                "-DCMAKE_RULE_MESSAGES=OFF",
                "--log-level=NOTICE",
                "-DENABLE_AVX=ON",
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

        bindir.mkdir(exist_ok=True)
        (Path(compiledir) / "build" / "strobealign").rename(binary)
    return binary, True


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
    binary, was_built = compile_strobealign(args.hash)


if __name__ == "__main__":
    main()
