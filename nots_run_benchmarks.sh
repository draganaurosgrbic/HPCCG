#!/bin/bash
#SBATCH --job-name=spmv_bench_final
#SBATCH --output=slurm_output_%j.log
#SBATCH --partition=commons
#SBATCH --constraint=sapphirerapids
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=192
#SBATCH --time=04:00:00

EXEC_PATH="$(pwd)/test_HPCCG"

THREADS=(1 2 4 8 12 16 24 32 36 48 64 96 128)
FORMATS=("csr" "ell7" "ell8" "tiled")
SIZES=("50 50 50" "100 100 100" "160 160 160")

export OMP_PLACES=cores

for sz in "${SIZES[@]}"; do
    echo "=========================================================="
    echo "MATRIX SIZE: $sz"
    echo "=========================================================="
    
    for fmt in "${FORMATS[@]}"; do
        for t in "${THREADS[@]}"; do
            
            if [ $t -ge 64 ]; then
                export OMP_PROC_BIND=spread
                RUN_CMD="numactl --interleave=all $EXEC_PATH"
            else
                export OMP_PROC_BIND=close
                RUN_CMD="$EXEC_PATH"
            fi
            
            export OMP_NUM_THREADS=$t
            
            echo "STARTING: Format=$fmt, Threads=$t, Size=$sz"
            $RUN_CMD $sz $fmt
            echo "----------------------------------------------------------"
        done
    done
done