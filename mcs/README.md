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

A pre-computed results file is available in this directory as `result-pre.csv`.
To reproduce the read alignment benchmark plots shown in the paper
(and the supplement), run:

    mkdir plots
    python plots.py -c config.yaml result-pre.csv plots/
    python plots.py --genome -c config.yaml result-pre.csv plots/


## Running the entire workflow

First, ensure the simulated datasets are available by running the workflow in
the root of this repository (`../Snakefile`). See the `README.md` file in
that directory. This will ensure that the symlinks to
`../datasets/` and `../genomes/` work.

Finally, run `snakemake -c 8`, where 8 is the number of threads to use. You can
use a different number, but our results were obtained with `-c 8`.
If you run this on a cluster, you should use
a Snakemake profile that specifies which resources to use, see
`../slurm/config.yaml` for an example.
With a profile, run `snakemake --profile=the-profile-dir/`
(where `the-profile-dir/` must contain a `config.yaml` file).


## Plots in the paper

The main output of the workflow is the CSV file `result.csv` and
the plots in `plots/`. The following plots are shown in the paper:

* `genome-CHM13-accuracy.pdf` -- Fig. 2 ("Accuracy of ... read alignment")
* `genome-CHM13-time.pdf` -- Fig. 3 ("Runtime of ... read alignment")
* `ends-se-accuracy.pdf` -- Suppl. Fig. 2 (Accuracy of single-end reads alignment)
* `ends-pe-accuracy.pdf` -- Suppl. Fig. 3 (Accuracy of paired-end reads alignment)
* `ends-se-time.pdf` -- Suppl. Fig. 4 (Runtime of single-end reads alignment)
* `ends-pe-time.pdf` -- Suppl. Fig. 5 (Runtime of paired-end reads alignment)
* `ends-se-aligned.pdf` -- Suppl. Fig. 6 (Mapping rate of single-end reads alignment)
* `ends-pe-aligned.pdf` -- Suppl. Fig. 7 (Mapping rate of paired-end reads alignment)
