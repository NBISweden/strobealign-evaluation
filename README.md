# Strobealign evaluation

This repository provides runnable workflows for evaluating
[strobealign](https://github.com/ksahlin/strobealign/).

The following workflows are available:

- `Snakefile`: Downloads reference genomes and generates libraries of simulated
  reads, grouped into three datasets named sim3, sim4, sim5 (with increasing
  mutaton rates). This workflow needs to be run as a prerequisite for the
  others.
- `minmax/Snakefile`: Reproduces the min/max evaluation experiments in the paper
  "Designing efficient randstrobes for sequence similarity
  analyses" by Karami et al., 2023. (see `minmax/README.md`)
- `simx/Snakefile`: Run strobealign, minimap2 and BWA-MEM on all datasets
  (sim3, sim4, sim5), measure accuracy, and plot the results.
- `mcs/Snakefile`: Reproduces the read alignment benchmark in the paper
  "Multi-context seeds enable fast and high-accuracy read mapping" by
  Tolstoganov et al., 2024 (see `mcs/README.md`).


## Datasets

The simulated datasets generated and used here are created in mostly
the same way as the datasets for the following paper:

Sahlin, K. Strobealign:
flexible seed size enables ultra-fast and accurate read alignment.
*Genome Biology 23, 260 (2022)*.
https://doi.org/10.1186/s13059-022-02831-7

See the Supplementary Information in that paper, section
"Note A: Simulations" for a description.

See also the original snakemake workflow file that documents the evaluation done
for the paper:
https://github.com/ksahlin/alignment_evaluation/blob/master/evaluation/Snakefile

Some differences:
- The paper mentions SIM3 and SIM4 datasets. We added a SIM5 dataset with high
  variability.
- The genomes in the paper are drosophila, maize, CHM13, rye. We added an
  E. coli pangenome (see below).
- The number of reads per library has been reduced to 1 million.


### *E. coli* pangenome

As fifth reference, an *E. coli* pangenome has been added
consisting of 50 randomly selected *E. coli* assemblies from RefSeq.
Because RefSeq changes, querying it is not reproducible.
Therefore, a pre-generated list of 100 *E. coli* assembly accessions
can be found in `ecoli-accessions.txt`.
It was generated in the following way:

    ncbi-genome-download --dry-run --assembly-levels complete --taxid 562 bacteria | sed 1d | cut -f1 | shuf | head -n 100 > ecoli-accessions.txt

The `ecoli50` datasets uses the first 50 accessions from that list.


### The number of reads has been reduced

The datasets in the Genome Biology paper contain 10 million simulated reads.
Here, we use only 1 million reads.

Note that there is a small difference between the results found by the workflow here
and the numbers reported in the Karami et al. paper
because two different ways were used to arrive at the smaller dataset:

- For the Karami et al. paper,
  the original 10 million simulated reads in each file were truncated to the
  first 1 million.
- The workflow in this repository simulates 1 million reads.

These datasets are not the same because the output generated by mason depends
on all input parameters,
and the number of reads it shall generate is one of these parameters.

To get the results to match exactly:

* Set `N_READS` in the `Snakefile` to 10'000'000
* Add a rule to truncate the generated data files to 1 million reads before
  mapping.

We have not made this modification as it adds hundreds of CPU hours for little gain.


# Running the workflow

Makes the required software available:

    conda env create --file environment.yml
    conda activate strobealign-eval

Download the genomes and simulate reads:

    snakemake --cores=all

This generates the folders `downloads/`, `genomes/` and `datasets/`.
The `downloads/` folder contains intermediate files and can be deleted.


## Issues

* `mason_variator` may fail with this error when running on maize or CHM13:

      /home/osboxes/seqan/apps/mason2/mason_variator.cpp:1161
      Assertion failed : snpRecord.getPos() != svRecord.getPos()
      should be true but was 0 (Should be generated non-overlapping (snp pos = 229421492, sv pos = 229421492).)

  We patch mason to work around this.
  An alternative is to compile mason without assertions.
