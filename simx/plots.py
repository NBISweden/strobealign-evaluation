#!/usr/bin/env python3
# If you run this script with `pipx run plots.py`, it will read the following
# inline script metadata and install the necessary dependencies into a
# temporary virtual enviroment
#
# /// script
# dependencies = [
#   "matplotlib",
#   "seaborn",
#   "pandas",
#   "pyyaml",
# ]
# ///

import os, sys
import argparse
import random
from pathlib import Path
from typing import Optional

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import yaml


MEASUREMENT_TYPES =  [
    # column name, label, logscale
    ("accuracy", "Accuracy (%)", False),
    ("aligned", "Percentage aligned", False),
    ("time", "Time (sec)", True),
    ("memory", "Memory usage (GB)", False),
]
LINEWIDTH = 3


def plot(
    table,
    palette,
    tools: dict[str, str],
    read_lengths,
    y: str,
    label: str,
    row: str = "genome",
    linewidth=1.5,
    xlim=(0, 500),
    logscale: bool = False,
    title: Optional[str] = None,
):
    g = sns.relplot(
        data=table,
        x="read_length",
        y=y,
        hue="tool",
        style="type",
        linewidth=linewidth,
        kind="line",
        col="dataset",
        row=row,
        facet_kws={"sharey": False},
        hue_order=tools,
        #col_order=["drosophila", "maize", "CHM13", "rye"],   # unused
        palette=palette,
        alpha=0.7,
    )
    g.figure.suptitle(title)
    g.set_axis_labels("Read length", label)
    legend_labels = ["Tool:"] + [name for key, name in tools.items()] + ["Type:", "align", "map"]
    sns.move_legend(
        g,
        loc="upper left",
        bbox_to_anchor=(0, 0),
        #facecolor="white",
        frameon=True,
        framealpha=1,
        ncols=len(legend_labels),
        labels=legend_labels,
    )

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
    return g


def configure(config_path):
    if config_path:
        with open(config_path) as f:
            config = yaml.safe_load(f)
    else:
        config = {}

    palette = {
        "minimap2": "tab:blue",
        "bwamem": "tab:orange",
        "xmapper": "red",
        # "strobealign_v071": "tab:green",
        # "strobealign_v0120_opt": "pink",
        # "strobealign_multicontext": "black",
    }

    for commit in config.get("commits", []):
        if "color" in commit:
            palette["strobealign-" + commit["key"]] = commit["color"]

    read_lengths = [50, 75, 100, 150, 200, 300, 500]
    # map short tool names to display names
    tools = {
        "minimap2": "minimap2",
        "bwamem": "BWA-MEM",
        "xmapper": "X-Mapper",
    }
    for commit in config["commits"]:
        tools["strobealign-" + commit["key"]] = commit["name"]

    return palette, read_lengths, tools


def plot_ends(df, outfolder, palette, read_lengths, tools, xlim, linewidth):
    for ends, table in df.groupby("ends"):
        title = "Single-end reads" if ends == "se" else "Paired-end reads"
        for y, label, logscale in MEASUREMENT_TYPES:
            fig = plot(
                table,
                palette,
                tools,
                read_lengths,
                y=y,
                logscale=logscale,
                row="genome",
                label=label,
                xlim=xlim,
                title=title,
                linewidth=linewidth,
            )
            if outfolder is not None:
                fig.savefig(outfolder / f"ends-{ends}-{y}.pdf")


def plot_genomes(df, outfolder, palette, read_lengths, tools, xlim, linewidth):
    for genome, table in df.groupby("genome"):
        for y, label, logscale in MEASUREMENT_TYPES:
            title = f"{label} – {genome}"
            fig = plot(
                table,
                palette,
                tools,
                read_lengths,
                y=y,
                logscale=logscale,
                row="ends",
                label=label,
                xlim=xlim,
                title=title,
                linewidth=linewidth,
            )
            if outfolder is not None:
                fig.savefig(outfolder / f"genome-{genome}-{y}.pdf")


def main(args):
    # Global plot settings
    #matplotlib.rcParams.update({"font.size": 18})
    sns.set(font_scale=1.4)
    sns.set_style("whitegrid")

    palette, read_lengths, tools = configure(args.config)
    xlim = (40, 510)

    table = pd.read_csv(args.csv)
    outfolder = Path(args.outfolder)
    if args.genome:
        plot_genomes(table, outfolder, palette, read_lengths, tools, xlim, LINEWIDTH)
    else:
        plot_ends(table, outfolder, palette, read_lengths, tools, xlim, LINEWIDTH)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--config", "-c", help="YAML configuration")
    parser.add_argument("--genome", action="store_true", help="Create genome-specific plots (with both single-end and paired-end measurements)")
    parser.add_argument("csv", help="Results file")
    parser.add_argument("outfolder", help="output folder")
    args = parser.parse_args()

    main(args)
