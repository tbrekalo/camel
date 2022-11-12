#! /usr/bin/fish
set -p fish_user_paths ./build/bin/
conda activate base

set threads 32
set work_dir /tmp/camel/
set ovlps /storage2/tbrekalo/yeast/diploid_sim/ovlps.paf
set reads /storage2/tbrekalo/yeast/diploid_sim/all_reads.fastq
set ref /storage2/tbrekalo/yeast/diploid_sim/NC_001133_9/ref.fa
set db_path ./eval.db

echo (which racon)
echo (which camel)

for e in (seq 0.1 0.1 0.1)
    for window_length in (seq 200 200 200)
        for exe in camel
            python3 ./scripts/pyeval --executable $exe \
                --args "$reads $ovlps --window-length $window_length --error-threshold $e" \
                --threads $threads \
                --work-dir $work_dir \
                --reference $ref \
                --database $db_path \
                --comment "grid search"
        end
    end
end
