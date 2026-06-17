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
    datasets: list[str],
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

    table = table[table["dataset"].isin(datasets)]
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
        facet_kws={"sharey": y == "time", "legend_out": False},
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
        sns.move_legend(
            g,
            loc="center left",
            bbox_to_anchor=(1, 0.5),
            labels=legend_labels,
            frameon=False,
        )
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

    # Fix a layout issue: If aligned percentage is very close to 100% for
    # all data points, the scaling of the y-axis is weird
    if y == "aligned":
        for ax in g.axes_dict.values():
            ylim = ax.get_ylim()
            if ylim[0] >= 99.9:
                ax.set_ylim((99.9, ylim[1]))

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
    }

    available_colors = ["tab:green", "black", "tab:red", "tab:purple", "tab:pink", "tab:brown", "tab:gray"]
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
            color = version["color"]
            try:
                available_colors.remove(version["color"])
            except ValueError:
                pass
        else:
            color = available_colors.pop(0)
        palette["strobealign-" + key] = color

    datasets = config["datasets"]
    modes = sorted(config.get("modes", ["align", "map"]))

    return palette, datasets, read_lengths, tools, modes


def plot_accuracy_table(table, tools, datasets, read_lengths, modes, title):
    ds_colors = ["#aeeaea", "#ffd6a0", "#b5e8b5", "#d4b8e0", "#f9c8c8"]

    filtered = table[table["dataset"].isin(datasets)]
    filtered = filtered[filtered["mode"].isin(modes)]
    filtered = filtered[filtered["read_length"].isin(read_lengths)]
    filtered = filtered[filtered["tool"].isin(tools.keys())]

    genomes = sorted(filtered["genome"].unique())
    if not genomes:
        return None

    tool_keys = list(tools.keys())
    tool_names = [tools[k] for k in tool_keys]
    all_col_tuples = [(ds, rl) for ds in datasets for rl in sorted(read_lengths)]

    cell_data = {
        (genome, ds, rl, tk): (
            filtered[(filtered["genome"] == genome) & (filtered["dataset"] == ds)
                     & (filtered["read_length"] == rl) & (filtered["tool"] == tk)]["accuracy"].mean()
            if not filtered[(filtered["genome"] == genome) & (filtered["dataset"] == ds)
                            & (filtered["read_length"] == rl) & (filtered["tool"] == tk)].empty
            else float("nan")
        )
        for genome in genomes for ds, rl in all_col_tuples for tk in tool_keys
    }

    col_tuples = [
        (ds, rl) for ds, rl in all_col_tuples
        if any(not pd.isna(cell_data[(g, ds, rl, tk)]) for g in genomes for tk in tool_keys)
    ]
    if not col_tuples:
        return None

    n_data_cols = len(col_tuples)
    rows_per_genome = 2 + len(tool_keys)
    fig_w = 1.2 + 2.8 + 0.62 * n_data_cols
    fig_h = 0.28 * len(genomes) * rows_per_genome + 0.5

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.axis("off")
    fig.suptitle(f"Accuracy (%) – {title}", fontsize=11, fontweight="bold")

    ds_color = {ds: ds_colors[i % len(ds_colors)] for i, ds in enumerate(datasets)}
    ds_first_col = {ds: min(ci for ci, (d, _) in enumerate(col_tuples) if d == ds) for ds in datasets if any(d == ds for d, _ in col_tuples)}

    for gi, genome in enumerate(genomes):
        block_top = fig_h - 0.35 - gi * rows_per_genome * 0.28
        block_bottom = block_top - rows_per_genome * 0.28

        best_per_col = {
            ci: max((v for tk in tool_keys if not pd.isna(v := cell_data[(genome, ds, rl, tk)])), default=None)
            for ci, (ds, rl) in enumerate(col_tuples)
        }

        def make_header_row(label_fn, bold_fn):
            text   = ["", ""] + [label_fn(ci, ds, rl) for ci, (ds, rl) in enumerate(col_tuples)]
            colors = ["#d0d0d0", "#e8e8e8"] + [ds_color[ds] for ds, rl in col_tuples]
            bold   = [False, False] + [bold_fn(ci, ds, rl) for ci, (ds, rl) in enumerate(col_tuples)]
            return text, colors, bold

        sim_text, sim_colors, sim_bold = make_header_row(
            lambda ci, ds, rl: ds if ds_first_col.get(ds) == ci else "",
            lambda ci, ds, rl: True,
        )
        rl_text, rl_colors, rl_bold = make_header_row(
            lambda ci, ds, rl: str(rl),
            lambda ci, ds, rl: False,
        )

        data_rows = []
        for ri, (tk, tname) in enumerate(zip(tool_keys, tool_names)):
            text   = ["" if ri > 0 else genome, tname]
            colors = ["#d0d0d0", "#e8e8e8"]
            bold   = [ri == 0, False]
            for ci, (ds, rl) in enumerate(col_tuples):
                v = cell_data[(genome, ds, rl, tk)]
                text.append(f"{v:.2f}" if not pd.isna(v) else "–")
                colors.append("white")
                bold.append(best_per_col[ci] is not None and not pd.isna(v) and abs(v - best_per_col[ci]) < 1e-9)
            data_rows.append((text, colors, bold))

        all_rows = [(sim_text, sim_colors, sim_bold), (rl_text, rl_colors, rl_bold)] + data_rows

        mpl_table = ax.table(
            cellText=[t for t, _, _ in all_rows],
            cellColours=[c for _, c, _ in all_rows],
            cellLoc="center",
            colWidths=[1.2 / fig_w, 2.8 / fig_w] + [0.62 / fig_w] * n_data_cols,
            bbox=[0.0, block_bottom / fig_h, 1.0, rows_per_genome * 0.28 / fig_h],
        )
        mpl_table.auto_set_font_size(False)
        mpl_table.set_fontsize(8)

        all_bold = [b for _, _, b in all_rows]
        for (row_idx, col_idx), cell in mpl_table.get_celld().items():
            cell.set_text_props(fontweight="bold" if all_bold[row_idx][col_idx] else "normal", fontsize=8)
            if col_idx <= 1:
                cell.get_text().set_ha("left")
            cell.set_linewidth(0.5)

    return fig


def plot_ends(df, outfolder, palette, datasets, read_lengths, tools, modes, xlogscale, linewidth, solid, legend):
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
                    datasets,
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
                    pdf.savefig(bbox_inches="tight")
                if outfolder is not None:
                    fig.savefig(outfolder / f"ends-{ends}-{y}.pdf")

            acc_fig = plot_accuracy_table(table, tools, datasets, read_lengths, modes, title)
            if acc_fig is not None:
                if outfolder is not None:
                    acc_fig.savefig(outfolder / f"ends-{ends}-accuracy-table.pdf", bbox_inches="tight")


def plot_genomes(df, outfolder, palette, datasets, read_lengths, tools, modes, xlogscale, linewidth, solid, legend):
    for genome, table in df.groupby("genome"):
        for y, label, logscale, _ in MEASUREMENT_TYPES:
            title = f"{label} – {genome}"
            fig = plot(
                table,
                palette,
                tools,
                modes,
                datasets,
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

    palette, datasets, read_lengths, tools, modes = configure(args.config)

    table = pd.read_csv(args.csv)
    table["time"] = table["time"] * 1E6 / table["read_count"]
    outfolder = Path(args.outfolder)
    if args.genome:
        plot_genomes(table, outfolder, palette, datasets, read_lengths, tools, modes, args.xlogscale, args.linewidth, args.solid, legend=args.legend)
    else:
        plot_ends(table, outfolder, palette, datasets, read_lengths, tools, modes, args.xlogscale, args.linewidth, args.solid, legend=args.legend)


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
