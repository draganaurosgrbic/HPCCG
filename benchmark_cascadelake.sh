#!/bin/bash
#SBATCH --job-name=spmv_cascadelake
#SBATCH --output=spmv_cascadelake_%j.log
#SBATCH --constraint=cascadelake
#SBATCH --partition=commons
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --time=04:00:00

THREADS=(1 2 4 8 12 16 20 24 32 40)
FORMATS=("csr" "ell7" "ell8" "tiled")
SIZES=("80 80 80" "100 100 100" "160 160 160" "200 200 200")

EXEC_PATH="$(pwd)/test_HPCCG"
export OMP_PLACES=cores

lscpu

for sz in "${SIZES[@]}"; do
    echo "=========================================================="
    echo "MATRIX SIZE: $sz"
    echo "=========================================================="
    
    for fmt in "${FORMATS[@]}"; do
        for t in "${THREADS[@]}"; do
            
            if [ $t -gt 20 ]; then
                export OMP_PROC_BIND=spread
                RUN_CMD="numactl --interleave=all $EXEC_PATH"
            else
                export OMP_PROC_BIND=close
                RUN_CMD="numactl --cpunodebind=0 --membind=0 $EXEC_PATH"
            fi
            
            export OMP_NUM_THREADS=$t
            echo "STARTING: Format=$fmt, Threads=$t, Size=$sz"
            $RUN_CMD $sz $fmt
            echo "----------------------------------------------------------"
        done
    done
done