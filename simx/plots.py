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
from contextlib import ExitStack

import matplotlib
from matplotlib.backends.backend_pdf import PdfPages

matplotlib.use("PDF")

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import yaml


MEASUREMENT_TYPES =  [
    # column name, label, logscale, include in single PDF
    ("accuracy", "Accuracy (%)", False, True),
    ("time", "Time (µs/read)", True, True),
#    ("saccuracy", "Score-based accuracy (%)", False),
    ("aligned", "Percentage aligned", False, False),
    ("memory", "Memory usage (GB)", False, False),
#    ("jaccuracy", "Jaccard accuracy (%)", False),
]


def plot(
    table,
    palette,
    tools: dict[str, str],
    modes: list[str],
    read_lengths,
    y: str,
    label: str,
    row: str = "genome",
    linewidth=1.5,
    solid: bool = False,
    logscale: bool = False,
    xlogscale: bool = True,
    title: Optional[str] = None,
    legend: str = "right",
):
    if y == "memory":
        tools = tools.copy()
        tools.pop("xmapper", None)
        table = table[table["tool"] != "xmapper"]

    table = table[table["mode"].isin(modes)]
    table = table[table["read_length"].isin(read_lengths)]

    use_style_for_tools = not solid and len(set(table["mode"])) == 1
    if use_style_for_tools and title is not None:
        mode = table["mode"].iloc[0]
        title += f" (mode: {mode})"
    g = sns.relplot(
        data=table,
        x="read_length",
        y=y,
        hue="tool",
        style="tool" if use_style_for_tools else "mode",
        linewidth=linewidth,
        kind="line",
        col="dataset",
        row=row,
        facet_kws={"sharey": y == "time"},
        hue_order=tools,
        #col_order=["drosophila", "maize", "CHM13", "rye"],   # unused
        palette=palette,
        alpha=0.7,
        clip_on=False,
        zorder=3,
    )
    g.figure.suptitle(title)
    g.set_axis_labels("Read length", label)

    if not use_style_for_tools:
        legend_labels = ["Tool:"] + [name for key, name in tools.items()]
        legend_labels += ["\nMode:"]
        assert modes == ["align"] or modes == ["map"] or modes == ["align", "map"]
        legend_labels += modes
    else:
        legend_labels = [name for key, name in tools.items()]

    if legend == "right":
        sns.move_legend(g, loc="right", labels=legend_labels)
    elif legend == "below":
        sns.move_legend(
            g,
            loc="upper left",
            bbox_to_anchor=(0, 0),
            frameon=True,
            ncols=len(legend_labels),
            labels=legend_labels,
        )

    if xlogscale:
        g.set(xscale="log")
    if logscale:
        g.set(yscale="log")
        yticks = []
        for m in [0.1, 1, 10, 100, 1000]:
            yticks.extend(m * n for n in range(1, 10))
        yticks.append(10000)
        low, high = table[y].min(), table[y].max()
        yticks = [t for t in yticks if low <= t <= high]
        g.set(yticks=yticks)

    read_lengths = sorted(set(read_lengths) & set(table["read_length"]))
    # Multiplying by 1.001 prevents 100 (etc.) from being displayed as
    # 10^2 in the tick labels
    xlim = (min(read_lengths), max(read_lengths) * 1.001)
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
        "xmapper": "turquoise",
        # "strobealign_v071": "tab:green",
        # "strobealign_v0120_opt": "pink",
        # "strobealign_multicontext": "black",
    }

    read_lengths = [50, 75, 100, 150, 200, 300, 500]
    if config["read-lengths"] is not None:
        read_lengths = sorted(config["read-lengths"])

    # map short tool names to display names
    names = {
        "minimap2": "minimap2",
        "bwamem": "BWA-MEM",
        "xmapper": "X-Mapper",
    }
    programs = config["programs"] if config["programs"] is not None else []
    tools = {program: names[program] for program in programs}

    for version in config.get("versions", []):
        key = version["commit"]
        if arg := version.get("arguments", ""):
            key += "-" + arg.replace(" ", "").replace("-", "").replace("=", "")
        tools["strobealign-" + key] = version["name"]
        if "color" in version:
            palette["strobealign-" + key] = version["color"]

    modes = sorted(config.get("modes", ["align", "map"]))

    return palette, read_lengths, tools, modes


def plot_ends(df, outfolder, palette, read_lengths, tools, modes, xlogscale, linewidth, solid, legend):
    with ExitStack() as stack:
        if outfolder is not None:
            pdf = stack.enter_context(PdfPages(outfolder / "ends.pdf"))
        else:
            pdf = None
        for ends, table in df.groupby("ends"):
            title = "Single-end reads" if ends == "se" else "Paired-end reads"
            for y, label, logscale, in_pdf in MEASUREMENT_TYPES:
                fig = plot(
                    table,
                    palette,
                    tools,
                    modes,
                    read_lengths,
                    y=y,
                    logscale=logscale,
                    xlogscale=xlogscale,
                    row="genome",
                    label=label,
                    title=title,
                    linewidth=linewidth,
                    solid=solid,
                    legend=legend,
                )
                if pdf is not None and in_pdf:
                    pdf.savefig()
                if outfolder is not None:
                    fig.savefig(outfolder / f"ends-{ends}-{y}.pdf")


def plot_genomes(df, outfolder, palette, read_lengths, tools, modes, xlogscale, linewidth, solid, legend):
    for genome, table in df.groupby("genome"):
        for y, label, logscale, _ in MEASUREMENT_TYPES:
            title = f"{label} – {genome}"
            fig = plot(
                table,
                palette,
                tools,
                modes,
                read_lengths,
                y=y,
                logscale=logscale,
                xlogscale=xlogscale,
                row="ends",
                label=label,
                title=title,
                linewidth=linewidth,
                solid=solid,
                legend=legend,
            )
            if outfolder is not None:
                fig.savefig(outfolder / f"genome-{genome}-{y}.pdf")


def main(args):
    # Global plot settings
    #matplotlib.rcParams.update({"font.size": 18})
    sns.set(font_scale=1.4)
    sns.set_style("whitegrid")

    palette, read_lengths, tools, modes = configure(args.config)

    table = pd.read_csv(args.csv)
    table["time"] = table["time"] * 1E6 / table["read_count"]
    outfolder = Path(args.outfolder)
    if args.genome:
        plot_genomes(table, outfolder, palette, read_lengths, tools, modes, args.xlogscale, args.linewidth, args.solid, legend=args.legend)
    else:
        plot_ends(table, outfolder, palette, read_lengths, tools, modes, args.xlogscale, args.linewidth, args.solid, legend=args.legend)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--config", "-c", help="YAML configuration")
    parser.add_argument("--linear-x", dest="xlogscale", default=True, action="store_false", help="Plot linear read length, not logscale")
    parser.add_argument("--linewidth", type=int, default=2)
    parser.add_argument("--solid", action="store_true", help="Line style depends on mode only (otherwise it depends on tool if plotting only one of map/align")
    parser.add_argument("--legend", choices=("below", "right"), default="right")
    parser.add_argument("--genome", action="store_true", help="Create genome-specific plots (with both single-end and paired-end measurements)")
    parser.add_argument("csv", help="Results file")
    parser.add_argument("outfolder", help="output folder")
    args = parser.parse_args()

    main(args)
