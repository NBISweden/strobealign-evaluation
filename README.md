# Strobealign min/max evaluation

This repository documents how to reproduce the min/max evaluation experiments
for the paper

"Designing efficient randstrobes for sequence similarity
analyses" by Karami et al., 2023.

(See "Implementing c_max in strobealign" in the appendix of that
paper.)


## Summary

To create a randstrobe, a syncmer is combined with a second one downstream.
This is done in a pseudorandom way by applying a function to the two syncmers
and choosing the combination that *minimizes* the output of that function.
An alternative is to *maximize* the function, but that gives similar results.

Here, we test what happens when running both the min and max version
and combining the results.


## Datasets

This evaluation uses the same drosophila, maize, CHM13 and rye simulated
datasets created for the following paper:

Sahlin, K. Strobealign:
flexible seed size enables ultra-fast and accurate read alignment.
*Genome Biol 23, 260 (2022)*.
https://doi.org/10.1186/s13059-022-02831-7

See the Supplementary Information in that paper, section
"Note A: Simulations" for a description of how to create the data.

Here, we only use the first 1 million reads of each simulated dataset.

The included snakemake workflow (see `Snakefile`) attempts to re-create the
simulated datasets (but see the notes below).


## Running the min/max evaluation


    conda env create --file environment.yml
    conda activate strobealign-eval
    snakemake -c 0

The final results are in `table.tex`.


## Issues

* `mason_variator` may fail with this error when running on maize or CHM13:

      /home/osboxes/seqan/apps/mason2/mason_variator.cpp:1161
      Assertion failed : snpRecord.getPos() != svRecord.getPos()
      should be true but was 0 (Should be generated non-overlapping (snp pos = 229421492, sv pos = 229421492).)

  This can be solved by using a self-compiled version of mason without assertions.

* The `Snakefile` included here exactly reproduces the steps for
  creating the simulated data, *except* that it only creates 1 million reads
  from the start instead of creating 10 million and taking the first 1 million.
  Since read simulation with mason is slow, this saves a lot of time,
  but because the set of simulated reads depends on the input parameters,
  the actual reads are different.


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
