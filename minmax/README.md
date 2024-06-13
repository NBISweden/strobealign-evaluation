# Min/max experiments

The workflow in this directory can be used to reproduce the min/max evaluation
experiments in the paper "Designing efficient randstrobes for sequence similarity
analyses" by Karami et al., 2023.

See also "Implementing c_max in strobealign" in the appendix of that paper.

For the paper,
[commit 5b26d727489b2b0c2177cac7a4f6e1193e1586bd](https://github.com/NBISweden/strobealign-evaluation/commit/5b26d727489b2b0c2177cac7a4f6e1193e1586bd).
was used. In that commit, the workflow consisted of a single `Snakefile`,
but that has now been split up into two separate workflows:
- `../Snakefile` (in the root of the repository) creates the datasets
- `Snakefile` (in this directory) runs the min/max evaluation as described
  below.

The two steps need to be run separately.

Only the `SIM3` dataset is used for the min/max evaluation.


# Summary

To create a randstrobe, a syncmer is combined with a second one downstream.
This is done in a pseudorandom way by applying a function to the two syncmers
and choosing the combination that *minimizes* the output of that function.
An alternative is to *maximize* the function, but that gives similar results.

Here, we test what happens when running both the min and max version
and combining the results.

The provided snakemake workflow (`Snakefile`) runs the evaluation "from scratch",
that is, no further input files need to be provided.

It does the following:
- Map reads with BWA-MEM
- Download and compile the two "min" and "max" strobealign versions
- Map reads with the two strobealign versions
- Combine min and max results
- Compute accuracy
- Create the result table in LaTeX format (see `table.tex`)

A pre-computed output file is provided in `precomputed-table.tex`.


## Running the min/max evaluation

Makes the required software available:

    conda env create --file environment.yml
    conda activate strobealign-eval

Download the genomes and simulate reads:

    snakemake --cores=all

Run the actual min/max evaluation:

    cd minmax
    snakemake --cores=all

The final results are in `minmax/table.tex`.
The file should be identical to the provided `minmax/precomputed-table.tex`.


## Commits

The strobealign commits used are:

- min: https://github.com/ksahlin/strobealign/commit/3223dc5946d9f38814e25a25149548dd146cc8d0 (original)
- max: https://github.com/ksahlin/strobealign/commit/6a837431f29fc3be3c3a74bea507538d9ea5abe7 (tag: min-max)


## Rules for combining reads

### Single-end reads

- If one read is unmapped, take the other
- If one read has lower alignment score, take the other

## Paired-end reads

- If one pair has lower sum of alignment scores, take the other
- If one pair is not marked as proper, take the other
- If one pair has fewer mapped reads, take the other
