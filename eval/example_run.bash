#!/bin/bash

exes=(camel racon)
threads=(64 32 16)
window_lengths=(200 500)

reads=(
  yeast_sim/diploid_sim/reads.fa
  yeast_S288C_simplex/reads.fastq
)

overlaps=(
  yeast_sim/diploid_sim/ovlps.paf
  yeast_S288C_simplex/ovlps.paf
)

references=(
  yeast_sim/diploid_sim/NC_001133_9/ref.fa
  yeast_S288C_simplex/GCF_000146045.2_R64_genomic.fna.gz
)

for i in $(seq 0 $((${#reads[@]} - 1)))
do
  for exe in "${exes[@]}"
  do
    job_desc=$(cat base_job.json | jq -c \
      --arg e $exe \
      --arg rp ${reads[$i]} \
      --arg op ${overlaps[$i]} \
      '.exe=$e | .readsPath=$rp | .overlapsPath=$op')

    echo $job_desc | jq ''
    python src/main.py <(echo $job_desc) -o evaluations/
  done
done
