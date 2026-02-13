#!/bin/bash
set -euo pipefail

# Run the pipeline for each configuration (YAML file)
# to create the .csv files
for config in chm13 xmapper strategies longreads hifi allgenomes memory; do
  echo "Configuration: ${config}"
  snakemake -p -c8 --configfile ${config}.yaml -- result.csv
  mv result.csv ${config}.csv
done
