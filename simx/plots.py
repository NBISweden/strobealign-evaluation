import os, sys
import argparse
import random
from typing import Optional

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


def plot(
    input_csv,
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
    indata = pd.read_csv(input_csv)
    g = sns.relplot(
        data=indata,
        x="read_length",
        y=y,
        hue="tool",
        style="type",
        linewidth=linewidth,
        kind="line",
        col="dataset",
        row="genome",
        facet_kws={"sharey": False},
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
    input_csv, outfolder, palette, tools, read_lengths, linewidth=2.5, xlim=(0, 500), title=None
):
    pdf_path = os.path.join(outfolder, "accuracy.pdf")
    plot(
        input_csv,
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
    input_csv, outfolder, palette, tools, read_lengths, linewidth=2.5, xlim=(0, 500), title=None
):
    pdf_path = os.path.join(outfolder, "aligned.pdf")
    plot(
        input_csv,
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
    input_csv, outfolder, palette, tools, read_lengths, linewidth=2.5, xlim=(0, 500), title=None
):
    pdf_path = os.path.join(outfolder, "memory.pdf")
    plot(
        input_csv,
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
    input_csv, outfolder, palette, tools, read_lengths, linewidth=2.5, xlim=(0, 500), title=None
):
    pdf_path = os.path.join(outfolder, "time.pdf")
    plot(
        input_csv,
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
    sns.set_style("whitegrid")
    palette = {
        "minimap2": "tab:blue",
        "bwamem": "tab:orange",
        "strobealign_v071": "tab:green",
        "strobealign_v0120_opt": "pink",
        "strobealign_multicontext": "black",
        "strobealign-main-3a97f6b": "tab:green",
        "strobealign-mcs-4ed851a": "black",
    }
    read_lengths = [50, 75, 100, 150, 200, 300, 500]
    tools = ["strobealign-main-3a97f6b", "strobealign-mcs-4ed851a", "bwamem"]
    xlim = (40, 260)
    csv = args.csv

    plot_accuracy(
        csv,
        args.outfolder,
        palette,
        tools,
        read_lengths,
        linewidth=2.5,
        xlim=xlim,
        title=args.title,
    )
    plot_percentage_aligned(
        csv,
        args.outfolder,
        palette,
        tools,
        read_lengths,
        linewidth=2.5,
        xlim=xlim,
        title=args.title,
    )
    plot_runtime(
        csv,
        args.outfolder,
        palette,
        tools,
        read_lengths,
        linewidth=2.5,
        xlim=xlim,
        title=args.title,
    )
    plot_memory_usage(
        csv,
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
        description="Calc identity",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--title", "-t", help="Optional title")
    parser.add_argument("csv", help="results file")
    parser.add_argument("outfolder", help="outfolder to plots.")
    args = parser.parse_args()

    main(args)
