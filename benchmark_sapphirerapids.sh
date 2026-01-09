#!/bin/bash
#SBATCH --job-name=spmv_sapphirerapids
#SBATCH --output=spmv_sapphirerapids_%j.log
#SBATCH --constraint=sapphirerapids
#SBATCH --partition=commons
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=96
#SBATCH --time=04:00:00

THREADS=(1 2 4 8 16 24 32 48 64 80 96)
FORMATS=("csr" "ell7" "ell8" "tiled")
SIZES=("100 100 100" "160 160 160" "250 250 250" "300 300 300")

EXEC_PATH="$(pwd)/test_HPCCG"
export OMP_PLACES=cores

lscpu

for sz in "${SIZES[@]}"; do
    echo "=========================================================="
    echo "MATRIX SIZE: $sz"
    echo "=========================================================="
    
    for fmt in "${FORMATS[@]}"; do
        for t in "${THREADS[@]}"; do
            
            if [ $t -gt 48 ]; then
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