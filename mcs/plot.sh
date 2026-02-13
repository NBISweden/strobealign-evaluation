#!/bin/bash
set -euo pipefail

function plot() {
  config="$1"
  shift

  nosim04=0
  if [[ $# -ge 1 && "$1" = "--nosim04" ]]; then
    nosim04=1
    shift
  fi
  echo "Plotting ${config} ..."
  if [ ${nosim04} = 1 ]; then
    actualcsv=${config}-nosim04.csv
    # Do not show noisy profile for SIM0 and SIM4
    grep -v 'k16s12l2u2,.*,sim[04]' ${config}.csv > ${actualcsv}
  else
    actualcsv=${config}.csv
  fi
  mkdir -p figures/${config}
  python plots.py "$@" -c ${config}.yaml ${actualcsv} figures/${config}/
  if [ ${nosim04} = 1 ]; then
    rm ${config}-nosim04.csv
  fi
}

#test -e figures && rm -r figures

# Figures 2 and 3
plot chm13 --nosim04 --genome
mv figures/chm13/genome-CHM13-{accuracy,time}.pdf figures/

# X-Mapper
plot xmapper --genome --solid
mv figures/xmapper/genome-CHM13-accuracy.pdf figures/xmapper-CHM13-accuracy.pdf
mv figures/xmapper/genome-CHM13-time.pdf figures/xmapper-CHM13-time.pdf

# Seed seed lookup strategies
plot strategies --genome --solid
mv figures/strategies/genome-CHM13-accuracy.pdf figures/strategies-CHM13-accuracy.pdf
mv figures/strategies/genome-CHM13-time.pdf figures/strategies-CHM13-time.pdf

# Long reads
plot longreads --solid --linear-x
mv figures/longreads/ends-se-time.pdf figures/longreads-time.pdf
mv figures/longreads/ends-se-accuracy.pdf figures/longreads-accuracy.pdf

# SIM1 and simulated HiFi reads
plot hifi --linear-x
mv figures/hifi/ends-se-accuracy.pdf figures/hifi-accuracy.pdf
mv figures/hifi/ends-se-time.pdf figures/hifi-runtime.pdf

# Full evaluation with all genomes
plot allgenomes --nosim04
mv figures/allgenomes/ends-se-accuracy.pdf figures/allgenomes-se-accuracy.pdf
mv figures/allgenomes/ends-pe-accuracy.pdf figures/allgenomes-pe-accuracy.pdf
mv figures/allgenomes/ends-se-time.pdf figures/allgenomes-se-time.pdf
mv figures/allgenomes/ends-pe-time.pdf figures/allgenomes-pe-time.pdf
mv figures/allgenomes/ends-se-aligned.pdf figures/allgenomes-se-aligned.pdf
mv figures/allgenomes/ends-pe-aligned.pdf figures/allgenomes-pe-aligned.pdf

plot memory
mv figures/memory/ends-pe-memory.pdf figures/memoryusage.pdf

echo "Done. Figures are in figures/"
