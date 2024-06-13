import os,sys
import argparse

import random

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except (ImportError, RuntimeError):
    print("COULD not import matplotlib")

import seaborn as sns
import pandas as pd
from matplotlib import pyplot


def plot_accuracy(input_csv, outfolder, palette, tools, read_lengths, linewidth = 2.5, xlim=(0,500)):
    matplotlib.rcParams.update({'font.size': 18})
    sns.set(font_scale=1.2)
    # sns.set(font_scale=1.6) # for main SIM3
    sns.set_style("whitegrid")
    indata = pd.read_csv(input_csv)
    # "minimap2",'accelalign', "bowtie2", "strobealign_map", "minimap2_map", "accelalign_map"
    # dashes = {"bwa_mem" : "", "bowtie2" : "",
    #             "minimap2" : "", "minimap2_map" : (5,5),
    #             "strobealign" : "", "strobealign_map" : (5,5),
    #             "accelalign" : "", "accelalign_map" : (5,5)}
    # print(indata)
    g = sns.relplot(data=indata, x="read_length", y="accuracy", hue="tool", style="type", linewidth = linewidth, kind="line", #dashes = dashes,
        col="dataset", hue_order = tools, facet_kws={'sharey': False}, # hue="datastructure", style="datastructure",
        col_wrap=2, col_order=["drosophila", "maize", "CHM13", "rye"], palette=palette)
         # col_order=["SIM3"], palette=palette)
    # ax = sns.lineplot(data=indata, x="k", y="unique", hue="datastructure", style="chr", palette = sns.color_palette()[:7])
    # axes = g.axes
    # g.set_titles("SIM3")
    g.set_axis_labels("Read length", "Accuracy")
    # g.set_xticklabels([18,24,30,36])
    # ax.set_ylabel("% unique")
    # ax.set_xlabel("k")
    # axes.set_xticks([18,24,30,36] )
    # ax.set_ylim((75, 100))
    # ylim=(80, 98),
    g.set( xlim=xlim, xticks=read_lengths)
    g.set_xticklabels(rotation=90, labels=read_lengths)
    g.tight_layout()
    # g.set(ylim=(95, 100))
    # ax.set_xticks([18,24,30,36])
    # plt.savefig(os.path.join(outfolder, "accuracy_plot.eps"))
    plt.savefig(os.path.join(outfolder, "accuracy_plot.pdf"))
    plt.close()

def plot_percentage_aligned(input_csv, outfolder, palette, tools, read_lengths, linewidth = 2.5, xlim=(0,500)):
    matplotlib.rcParams.update({'font.size': 18})
    sns.set(font_scale=1.2)
    # sns.set(font_scale=1.6) # for main SIM3
    sns.set_style("whitegrid")

    indata = pd.read_csv(input_csv)

    # dashes = { "strobemap" : "", 
    #             "spaced_sparse": (5,5),
    #             "spaced_dense": (5,5),
    #              "minstrobes2" : (1,1),
    #              "minstrobes3" : (1,1),
    #              "randstrobes2" : (1,1),
    #              "randstrobes3" : (1,1),
    #              "hybridstrobes2" : (1,1),
    #              "hybridstrobes3" : (1,1)}
    # print(indata)
    g = sns.relplot(
        data=indata, x="read_length", y="aligned", hue="tool", style="type", linewidth = linewidth,
        col="dataset", kind="line",  hue_order = tools, #dashes = dashes, hue="datastructure", style="datastructure",
        col_wrap=2, col_order=["drosophila", "maize", "CHM13", "rye"], palette=palette)
         # col_order=["SIM3"], palette=palette)
    # ax = sns.lineplot(data=indata, x="k", y="unique", hue="datastructure", style="chr", palette = sns.color_palette()[:7])
    # axes = g.axes
    # g.set_titles("SIM3")
    g.set_axis_labels("Read length", "Percentage aligned")
    # g.set_xticklabels([18,24,30,36])
    # ax.set_ylabel("% unique")
    # ax.set_xlabel("k")
    # axes.set_xticks([18,24,30,36] )
    # ax.set_ylim((75, 100))
    g.set(ylim=(95, 100), xlim=xlim, xticks=read_lengths)
    g.set_xticklabels(rotation=90, labels=read_lengths)
    g.tight_layout()
    # g.set(ylim=(95, 100))
    # ax.set_xticks([18,24,30,36])
    # plt.savefig(os.path.join(outfolder, "percentage_aligned_plot.eps"))
    plt.savefig(os.path.join(outfolder, "percentage_aligned_plot.pdf"))
    plt.close()

def plot_overaligned(input_csv, outfolder, palette, tools, read_lengths, linewidth = 2.5, xlim=(0,500)):
    matplotlib.rcParams.update({'font.size': 18})
    sns.set(font_scale=1.2)
    # sns.set(font_scale=1.6) # for main SIM3
    sns.set_style("whitegrid")

    indata = pd.read_csv(input_csv)

    # dashes = { "strobemap" : "", 
    #             "spaced_sparse": (5,5),
    #             "spaced_dense": (5,5),
    #              "minstrobes2" : (1,1),
    #              "minstrobes3" : (1,1),
    #              "randstrobes2" : (1,1),
    #              "randstrobes3" : (1,1),
    #              "hybridstrobes2" : (1,1),
    #              "hybridstrobes3" : (1,1)}
    # print(indata)
    g = sns.relplot(
        data=indata, x="read_length", y="overaligned", hue="tool", style="type", linewidth = linewidth,
        col="dataset", kind="line", hue_order = tools,  #dashes = dashes, hue="datastructure", style="datastructure",
        col_wrap=2, col_order=["drosophila", "maize", "CHM13", "rye"], palette=palette)
         # col_order=["SIM3"], palette=palette)
    # ax = sns.lineplot(data=indata, x="k", y="unique", hue="datastructure", style="chr", palette = sns.color_palette()[:7])
    # axes = g.axes
    # g.set_titles("SIM3")
    g.set_axis_labels("Read length", "Overaligned")
    # g.set_xticklabels([18,24,30,36])
    # ax.set_ylabel("% unique")
    # ax.set_xlabel("k")
    # axes.set_xticks([18,24,30,36] )
    # ax.set_ylim((75, 100))
    g.set(xlim=xlim, xticks=read_lengths) #ylim=(40, 100),
    g.set_xticklabels(rotation=90, labels=read_lengths)
    g.tight_layout()
    # g.set(ylim=(95, 100))
    # ax.set_xticks([18,24,30,36])
    # plt.savefig(os.path.join(outfolder, "overaligned_plot.eps"))
    plt.savefig(os.path.join(outfolder, "overaligned_plot.pdf"))
    plt.close()


def plot_memory_usage(input_csv, outfolder, palette, tools, read_lengths, linewidth = 2.5, xlim=(0,500)):
    matplotlib.rcParams.update({'font.size': 18})
    sns.set(font_scale=1.2)
    # sns.set(font_scale=1.6) # for main SIM3
    # tool,dataset,read_length,time,memory
    indata = pd.read_csv(input_csv)

    sns.set_style("whitegrid")

    g = sns.relplot(
        data=indata, x="read_length", y="memory", hue="tool", style="type", linewidth = linewidth, 
        col="dataset", kind="line", hue_order = tools,  #dashes = dashes, hue="datastructure", style="datastructure",
        col_wrap=2, col_order=["drosophila", "maize", "CHM13", "rye"], palette=palette)
         # col_order=["SIM3"], palette=palette)
    # ax = sns.lineplot(data=indata, x="k", y="unique", hue="datastructure", style="chr", palette = sns.color_palette()[:7])
    # axes = g.axes
    # g.set_titles("SIM3")
    g.set_axis_labels("Read length", "Memory usage (Gb)")
    # g.set_xticklabels([18,24,30,36])
    # ax.set_ylabel("% unique")
    # ax.set_xlabel("k")
    # axes.set_xticks([18,24,30,36] )
    # ax.set_ylim((75, 100))
    g.set(xlim=xlim, xticks=read_lengths) #ylim=(40, 100),
    g.set_xticklabels(rotation=90, labels=read_lengths)
    g.tight_layout()
    # g.set(ylim=(95, 100))
    # ax.set_xticks([18,24,30,36])
    # plt.savefig(os.path.join(outfolder, "memory_plot.eps"))
    plt.savefig(os.path.join(outfolder, "memory_plot.pdf"))
    plt.close()

def plot_runtime(input_csv, outfolder, palette, tools, read_lengths, linewidth = 2.5, xlim=(0,500)):
    matplotlib.rcParams.update({'font.size': 18})
    sns.set(font_scale=1.2)
    # sns.set(font_scale=1.6) # for main SIM3
    # tool,dataset,read_length,time,memory
    sns.set_style("whitegrid")

    indata = pd.read_csv(input_csv)
    g = sns.relplot(
        data=indata, x="read_length", y="time", hue="tool", style="type", linewidth = linewidth,
        col="dataset", kind="line", hue_order = tools,  #dashes = dashes, hue="datastructure", style="datastructure",
        col_wrap=2, col_order=["drosophila", "maize", "CHM13", "rye"], palette=palette)
         # col_order=["SIM3"], palette=palette) # for main SIM3
    # ax = sns.lineplot(data=indata, x="k", y="unique", hue="datastructure", style="chr", palette = sns.color_palette()[:7])
    # axes = g.axes
    # g.set_titles("SIM3")
    g.set_axis_labels("Read length", "Time (sec)")
    # g.set_xticklabels([18,24,30,36])
    # ax.set_ylabel("% unique")
    # ax.set_xlabel("k")
    # axes.set_xticks([18,24,30,36] )
    # ax.set_ylim((75, 100))
    g.set(yscale="log")
    g.set( yticks= [i for i in range(10,99,10)] + [i for i in range(100,999,100)] + [i for i in range(1000,1999,1000)]) # + [i for i in range(10000,39999,10000)]) #, ylim=(0, 5200))
    # g.set_yticklabels( ["100"] + ["" for i in range(200,999,100)] + ["1000"] +  ["" for i in range(2000,9999,1000)] + ["10000"] +  [i for i in range(1000,2999,1000)]])

    g.set(xlim=xlim, xticks=read_lengths) #ylim=(40, 100),
    g.set_xticklabels(rotation=90, labels=read_lengths)
    g.tight_layout()
    # g.set( yticks=[0,1000,2000,4000,6000,12000,18000,24000]) #ylim=(40, 100),

    # g.set(ylim=(95, 100))
    # ax.set_xticks([18,24,30,36])
    # plt.savefig(os.path.join(outfolder, "time_plot.eps"))
    plt.savefig(os.path.join(outfolder, "time_plot.pdf"))
    plt.close()

def add_column(infile):
    mod_outfile = open(infile + "_mod.csv", "w")
    for i, line in enumerate(open(infile,'r')):
        if i == 0:
            line = line.strip()+",type\n"
            mod_outfile.write(line)
            continue
        vals = line.strip().split(",")
        is_aln = True
        if "_map" in vals[0]:
            v_tmp = vals[0][:-4]
            vals[0] = v_tmp
            is_aln = False

        if is_aln:
            vals.append("align")
        else:
            vals.append("map")

        mod_line = ",".join([v for v in vals]) + "\n"
        mod_outfile.write(mod_line)
    mod_outfile.close()
    return mod_outfile.name



def main(args):
    sns.set_style("whitegrid")
    # palette = {
    # 'minimap2': 'tab:blue',
    # 'strobealign_v071': 'tab:green',
    # 'bwa_mem': 'tab:orange',
    # 'accelalign': 'tab:red',
    # 'bowtie2' : 'tab:purple',
    # 'urmap' : 'tab:grey',
    # 'snap' : 'tab:pink',
    # "bwa_mem2" : 'black',
    # "strobealign_v080" : 'cyan',
    # "strobealign_master_preindexed" : 'gold'
    # }
    # tools =["minimap2", "bwa_mem", 'accelalign', "bowtie2", "snap", "bwa_mem2", "strobealign_v071", "strobealign_master_preindexed", "strobealign_v080"] #, "strobealign_mixed"] "urmap",

    palette = {
    'minimap2': 'tab:blue',
    'bwa_mem': 'tab:orange',
    'strobealign_v071': 'tab:green',
    "strobealign_v0120_opt" : 'pink',
    "strobealign_multicontext" : 'black'
    }
    # read_lengths = [50, 75, 91, 100, 111, 125, 136, 150, 176, 200, 250, 300, 500]
    read_lengths = [50, 75, 100, 150, 200, 250]
    tools = ["minimap2", "bwa_mem", "strobealign_v071", "strobealign_v0120_opt", "strobealign_multicontext"]
    xlim=(40,260)
    accuracy_csv = add_column(args.accuracy_csv)
    runtime_mem_csv = add_column(args.runtime_mem_csv)

    plot_accuracy(accuracy_csv, args.outfolder, palette, tools, read_lengths, linewidth = 2.5, xlim=xlim)
    plot_percentage_aligned(accuracy_csv, args.outfolder, palette, tools, read_lengths, linewidth = 2.5, xlim=xlim)
    # plot_overaligned(accuracy_csv, args.outfolder, palette, tools, read_lengths, linewidth = 2.5, xlim=xlim)
    plot_runtime(runtime_mem_csv, args.outfolder, palette, tools, read_lengths, linewidth = 2.5, xlim=xlim)
    plot_memory_usage(runtime_mem_csv, args.outfolder, palette, tools, read_lengths, linewidth = 2.5, xlim=xlim)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calc identity", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('accuracy_csv', type=str, help='results file')
    parser.add_argument('runtime_mem_csv', type=str, help='results file')
    parser.add_argument('outfolder', type=str,  help='outfolder to plots.')
    args = parser.parse_args()



    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()

    main(args)