# Read alignment benchmark for multi-context seeds

The Snakemake workflow in this directory reproduces the read alignment
benchmark in the paper

"Multi-context seeds enable fast and high-accuracy read mapping"
by Tolstoganov et al.

To get the workflow exactly as it was run for the pre-print at
<https://doi.org/10.1101/2024.10.29.620855>, please use commit
569827b28b76b643cb975a17e001d426067264cd. Some changes have been
made to the workflow since then.


## Reproducing only the plots

Pre-computed results files in CSV format are available in this directory.
To reproduce the read alignment benchmark plots shown in the paper
(main text and supplement), run `./plot.sh` in an activated Conda environment
(see `README.md` in the root directory for how to create the environment.)
Plots are placed into `figures/`.

## Running the entire workflow

First, ensure the simulated datasets are available by running the workflow in
the root of this repository (`../Snakefile`). See the `README.md` file in
that directory. This will ensure that the symlinks to
`../datasets/` and `../genomes/` work.

Next, run the shell script `./run.sh`. This runs snakemake once on
each of the various benchmark configurations (`.yaml` files in this directory)
in order to create a CSV file for each one.

The script runs snakemake using eight cores (`-c 8`). We ran our benchmarks on
a dedicated PC with very few other processes running in order to get reliable
runtime measurements. The machine has 8 cores and 64 GB RAM. However, it will
take a long time to run everything from scratch (multiple days).

To run the workflow on a cluster, the `./run.sh` script should be simple to
modify. If you use SLURM, you can use the Snakemake profile in
`../slurm/config.yaml` as a template. (Use `snakemake --profile=the-profile-dir/`
to run snakemake with a profile.)

After running the workflow, create the final plots by running `./plot.sh` as
described in the previous section.


## Plots in the paper

The main output of the workflow are the CSV files in this directory and
the plots in `figures/`. Plots in subdirectories of `figures/` (such as
`figures/longreads/`) are created by the workflow, but not shown in the
paper.
