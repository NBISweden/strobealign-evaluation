#!/bin/bash
set -euo pipefail

function run_and_plot() {
  config="$1"
  shift
  echo "Configuration: ${config}"
  snakemake -p -c8 --configfile ${config}.yaml -- result.csv
  mv result.csv ${config}.csv
  mkdir -p plots-${config}
  echo Plotting ...
  python plots.py "$@" -c ${config}.yaml ${config}.csv plots-${config}/
}

mkdir -p figures

# Figures 2 and 3
run_and_plot chm13 --legend=below --genome
mv plots-chm13/genome-CHM13-{accuracy,time}.pdf figures/

# X-Mapper
run_and_plot xmapper --genome --solid
mv plots-xmapper/genome-CHM13-accuracy.pdf figures/xmapper-CHM13-accuracy.pdf
mv plots-xmapper/genome-CHM13-time.pdf figures/xmapper-CHM13-time.pdf

# Seed selection strategies
run_and_plot strategies
run_and_plot strategies --genome --solid
mv plots-strategies/genome-CHM13-accuracy.pdf figures/strategies-CHM13-accuracy.pdf
mv plots-strategies/genome-CHM13-time.pdf figures/strategies-CHM13-time.pdf
mv plots-strategies/ends-se-time.pdf figures/strategies-se-time.pdf
mv plots-strategies/ends-pe-time.pdf figures/strategies-pe-time.pdf
mv plots-strategies/ends-se-accuracy.pdf figures/strategies-se-accuracy.pdf
mv plots-strategies/ends-pe-accuracy.pdf figures/strategies-pe-accuracy.pdf

# Long reads
run_and_plot longreads --solid --linear-x
mv plots-longreads/ends-se-time.pdf figures/longreads-time.pdf
mv plots-longreads/ends-se-accuracy.pdf figures/longreads-accuracy.pdf

# Full evaluation with all genomes ("data dump")
run_and_plot allgenomes
mv plots-allgenomes/ends-se-accuracy.pdf figures/allgenomes-se-accuracy.pdf
mv plots-allgenomes/ends-pe-accuracy.pdf figures/allgenomes-pe-accuracy.pdf
mv plots-allgenomes/ends-se-time.pdf figures/allgenomes-se-time.pdf
mv plots-allgenomes/ends-pe-time.pdf figures/allgenomes-pe-time.pdf
mv plots-allgenomes/ends-se-aligned.pdf figures/allgenomes-se-aligned.pdf
mv plots-allgenomes/ends-pe-aligned.pdf figures/allgenomes-pe-aligned.pdf

run_and_plot memory
mv plots-memory/ends-pe-memory.pdf figures/memoryusage.pdf
