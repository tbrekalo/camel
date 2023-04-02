# evaluation

This folder contains scripts for evaluating read correction tools, racon and camel. `src/main.py` is the main entry point.

```bash
usage: evaluation [-h] -r REFERENCE -o OUTPUT [-t THREADS] config

evaluate read correction tool

positional arguments:
  config                utf-8 json file with tool configuration

options:
  -h, --help            show this help message and exit
  -r REFERENCE, --reference REFERENCE
                        reference reads
  -o OUTPUT, --output OUTPUT
                        output folder containing runtime information
  -t THREADS, --threads THREADS
                        nmber of threads used for evaluation
```

`base_job.json` is a template for creating json job descriptions. `example_run.bash` illustrates iterating over sample data with different job configurations.

```bash
#!/bin/bash

exes=(camel racon)
threads=(64 32 16)
window_lengths=(200 500)
base_out_dir=evaluations/

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

      run_name=$exe\_$(date +"%Y-%m-%d_%H-%M-%S")
      run_dir=$base_out_dir/$run_name

      echo $run_dir
      mkdir $run_dir

    echo $job_desc | jq ''
    python src/main.py <(echo $job_desc) \
      -r ${references[$i]} \
      -o $run_dir \
      -t 32
  done
done
```

Each run generates a directory with benchmarks along with corrected sequences.

```bash
eval_out/
└── camel_23-03-28_04-20
    ├── info.json
    └── reads.fa
```

info.json

```json
{
  "taskConfig": {
    "exe": "camel",
    "threads": 32,
    "windowLength": 200,
    "errorThreshold": 0.3,
    "readsPath": "/home/tbrekalo/camel/data/example.fq",
    "overlapsPath": "/home/tbrekalo/camel/data/example.paf"
  },
  "taskRun": {
    "peakMemoryMib": 43.03,
    "runtimeS": 0.18,
    "accuracies": [0.99, 0.97, ...]
  }
}
```
