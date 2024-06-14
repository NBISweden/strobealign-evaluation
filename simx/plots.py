import os, sys
import argparse
import random

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
        kind="line",  # dashes = dashes,
        col="dataset",
        row="genome",
        ###hue_order=tools,  ### TODO
        facet_kws={"sharey": False},  # hue="datastructure", style="datastructure",
        # col_wrap=2,
        ###col_order=["drosophila", "maize", "CHM13", "rye"],   ### TODO
        ###palette=palette,   ### TODO
    )
    g.set_axis_labels("Read length", label)
    if logscale:
        g.set(yscale="log")
        g.set(
            yticks=[i for i in range(10, 99, 10)]
            + [i for i in range(100, 999, 100)]
            + [i for i in range(1000, 1999, 1000)]
        )  # + [i for i in range(10000,39999,10000)]) #, ylim=(0, 5200))

    g.set(xlim=xlim, xticks=read_lengths)
    g.set_xticklabels(rotation=90, labels=read_lengths)
    g.tight_layout()
    plt.savefig(pdf_path)
    plt.close()


def plot_accuracy(
    input_csv, outfolder, palette, tools, read_lengths, linewidth=2.5, xlim=(0, 500)
):
    pdf_path = os.path.join(outfolder, "accuracy_plot.pdf")
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
    )


def plot_percentage_aligned(
    input_csv, outfolder, palette, tools, read_lengths, linewidth=2.5, xlim=(0, 500)
):
    pdf_path = os.path.join(outfolder, "percentage_aligned_plot.pdf")
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
    )


def plot_memory_usage(
    input_csv, outfolder, palette, tools, read_lengths, linewidth=2.5, xlim=(0, 500)
):
    pdf_path = os.path.join(outfolder, "memory_plot.pdf")
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
    )


def plot_runtime(
    input_csv, outfolder, palette, tools, read_lengths, linewidth=2.5, xlim=(0, 500)
):
    pdf_path = os.path.join(outfolder, "time_plot.pdf")
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
    )


def add_column(infile):
    with open(infile + "_mod.csv", "w") as mod_outfile:
        for i, line in enumerate(open(infile, "r")):
            if i == 0:
                line = line.strip() + ",type\n"
                mod_outfile.write(line)
                continue
            vals = line.strip().split(",")
            is_aln = True
            if "_map" in vals[0]:
                v_tmp = vals[0][:-4]
                vals[0] = v_tmp
                is_aln = False

            vals.append("align" if is_aln else "map")
            mod_line = ",".join(vals) + "\n"
            mod_outfile.write(mod_line)

    return mod_outfile.name


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
    csv = add_column(args.csv)

    plot_accuracy(
        csv,
        args.outfolder,
        palette,
        tools,
        read_lengths,
        linewidth=2.5,
        xlim=xlim,
    )
    plot_percentage_aligned(
        csv,
        args.outfolder,
        palette,
        tools,
        read_lengths,
        linewidth=2.5,
        xlim=xlim,
    )
    plot_runtime(
        csv,
        args.outfolder,
        palette,
        tools,
        read_lengths,
        linewidth=2.5,
        xlim=xlim,
    )
    plot_memory_usage(
        csv,
        args.outfolder,
        palette,
        tools,
        read_lengths,
        linewidth=2.5,
        xlim=xlim,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calc identity",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("csv", help="results file")
    parser.add_argument("outfolder", help="outfolder to plots.")
    args = parser.parse_args()

    main(args)
