# Development utilities

## Evaluation

`eval` folder contains python evaluator and ORM modeling for string data used in report construction in `report.ipynb`. Below is an example fish script for running the evaluator in multiple parameters storing the input parameters along with evaluation results in sqlite3 database file. ``

```bash
#! /usr/bin/fish
set -p fish_user_paths ./build/bin/
conda activate base # puts racon executable in path

set work_dir /tmp/camel/
set S288C_Nanopore_R7.fastq.gz
set ovlps.paf
set ref reference.fasta.gz
set db_path ./eval.db

echo (which racon)
echo (which camel)

for threads in (seq 32 16 64)
    for window_length in 100 200 500 1000 2000
        for exe in camel racon
            python3 ./eval/evaluator --executable $exe \
                --args "$reads $ovlps --window-length $window_length" \
                --threads $threads \
                --work-dir $work_dir \
                --reference $ref \
                --database $db_path \
                --comment "benchmark"
        end
    end
end

