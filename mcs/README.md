# Read alignment benchmark for multi-context seeds

The Snakemake workflow in this directory reproduces the read alignment
benchmark in the paper

"Multi-context seeds enable fast and high-accuracy read mapping"
by Tolstoganov et al., 2024.

## Reproducing only the plots

A pre-computed results file is available in this directory as `result.csv`.
To reproduce the read alignment benchmark plots shown in the paper
(and the supplement), run:

    mkdir plots
    python plots.py -c config.yaml result.csv plots/
    python plots.py --genome -c config.yaml result.csv plots/


## Running the entire workflow

First, ensure the simulated datasets are available by running the workflow in
the root of this repository (`../Snakefile`). See the `README.md` file in
that directory. This will ensure that the symlinks to
`../datasets/` and `../genomes/` work.

Finally, run `snakemake -c all`. If you run this on a cluster, you should use
a Snakemake profile that specifies which resources to use, see
`../slurm/config.yaml` for an example.
With a profile, run `snakemake --profile=the-profile-dir/`
(where `the-profile-dir/` must contain a `config.yaml` file).


## Plots in the paper

The main output of the workflow is the CSV file `results/result.csv` and
the plots in `plots/`. The following plots are shown in the paper:

* `genome-CHM13-accuracy.pdf` -- Fig. 2 ("Accuracy of ... read alignment")
* `genome-CHM13-time.pdf` -- Fig. 3 ("Runtime of ... read alignment")
* `ends-se-accuracy.pdf` -- Suppl. Fig. 2 (Accuracy of single-end reads alignment)
* `ends-pe-accuracy.pdf` -- Suppl. Fig. 3 (Accuracy of paired-end reads alignment)
* `ends-se-time.pdf` -- Suppl. Fig. 4 (Runtime of single-end reads alignment)
* `ends-pe-time.pdf` -- Suppl. Fig. 5 (Runtime of paired-end reads alignment)
* `ends-se-aligned.pdf` -- Suppl. Fig. 6 (Mapping rate of single-end reads alignment)
* `ends-pe-aligned.pdf` -- Suppl. Fig. 7 (Mapping rate of paired-end reads alignment)
