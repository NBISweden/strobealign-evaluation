import os, sys
import argparse
import random
from typing import Optional

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import yaml


def plot(
    table,
    pdf_path,
    palette,
    tools,
    read_lengths,
    y: str,
    label: str,
    linewidth=2.5,
    xlim=(0, 500),
    logscale: bool = False,
    title: Optional[str] = None,
):
    matplotlib.rcParams.update({"font.size": 18})
    sns.set(font_scale=1.2)
    sns.set_style("whitegrid")
    g = sns.relplot(
        data=table,
        x="read_length",
        y=y,
        hue="tool",
        style="type",
        linewidth=linewidth,
        kind="line",
        col="dataset",
        row="genome",
        facet_kws={"sharey": False},
        hue_order=tools,
        #col_order=["drosophila", "maize", "CHM13", "rye"],   # TODO
        palette=palette,
    )
    g.figure.suptitle(title)
    g.set_axis_labels("Read length", label)
    if logscale:
        g.set(yscale="log")
        g.set(
            yticks=list(range(10, 99, 10))
            + list(range(100, 999, 100))
            + list(range(1000, 1999, 1000))
        )  # + list(range(10000,39999,10000))

    g.set(xlim=xlim, xticks=read_lengths)
    g.set_xticklabels(rotation=90, labels=read_lengths)
    g.tight_layout()
    plt.savefig(pdf_path)
    plt.close()


def plot_accuracy(
    table, outfolder, palette, tools, read_lengths, linewidth=2.5, xlim=(0, 500), title=None
):
    pdf_path = os.path.join(outfolder, "accuracy.pdf")
    plot(
        table,
        pdf_path,
        palette,
        tools,
        read_lengths,
        y="accuracy",
        label="Accuracy (%)",
        linewidth=linewidth,
        xlim=xlim,
        title=title,
    )


def plot_percentage_aligned(
    table, outfolder, palette, tools, read_lengths, linewidth=2.5, xlim=(0, 500), title=None
):
    pdf_path = os.path.join(outfolder, "aligned.pdf")
    plot(
        table,
        pdf_path,
        palette,
        tools,
        read_lengths,
        y="aligned",
        label="Percentage aligned",
        linewidth=linewidth,
        xlim=xlim,
        title=title,
    )


def plot_memory_usage(
    table, outfolder, palette, tools, read_lengths, linewidth=2.5, xlim=(0, 500), title=None
):
    pdf_path = os.path.join(outfolder, "memory.pdf")
    plot(
        table,
        pdf_path,
        palette,
        tools,
        read_lengths,
        y="memory",
        label="Memory usage (GB)",
        linewidth=linewidth,
        xlim=xlim,
        title=title,
    )


def plot_runtime(
    table, outfolder, palette, tools, read_lengths, linewidth=2.5, xlim=(0, 500), title=None
):
    pdf_path = os.path.join(outfolder, "time.pdf")
    plot(
        table,
        pdf_path,
        palette,
        tools,
        read_lengths,
        y="time",
        label="Time (sec)",
        linewidth=linewidth,
        xlim=xlim,
        logscale=True,
        title=title,
    )


def main(args):
    if args.config:
        with open(args.config) as f:
            config = yaml.safe_load(f)
    else:
        config = {}

    sns.set_style("whitegrid")
    palette = {
        "minimap2": "tab:blue",
        "bwamem": "tab:orange",
        # "strobealign_v071": "tab:green",
        # "strobealign_v0120_opt": "pink",
        # "strobealign_multicontext": "black",
    }

    for commit in config.get("commits", []):
        if "color" in commit:
            palette["strobealign-" + commit["name"]] = commit["color"]

    read_lengths = [50, 75, 100, 150, 200, 300, 500]
    tools = [
        "minimap2",
        "bwamem",
    ]
    for commit in config["commits"]:
        tools.append("strobealign-" + commit["name"])

    xlim = (40, 260)
    csv = args.csv


def main(args):
    table = pd.read_csv(args.csv)
    plot_accuracy(
        table,
        args.outfolder,
        palette,
        tools,
        read_lengths,
        linewidth=2.5,
        xlim=xlim,
        title=args.title,
    )
    plot_percentage_aligned(
        table,
        args.outfolder,
        palette,
        tools,
        read_lengths,
        linewidth=2.5,
        xlim=xlim,
        title=args.title,
    )
    plot_runtime(
        table,
        args.outfolder,
        palette,
        tools,
        read_lengths,
        linewidth=2.5,
        xlim=xlim,
        title=args.title,
    )
    plot_memory_usage(
        table,
        args.outfolder,
        palette,
        tools,
        read_lengths,
        linewidth=2.5,
        xlim=xlim,
        title=args.title,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--title", "-t", help="Optional title")
    parser.add_argument("--config", "-c", help="YAML configuration")
    parser.add_argument("csv", help="results file")
    parser.add_argument("outfolder", help="output folder")
    args = parser.parse_args()

    main(args)
