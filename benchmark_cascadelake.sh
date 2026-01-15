#!/bin/bash
#SBATCH --job-name=spmv_cascadelake
#SBATCH --output=spmv_cascadelake_%j.log
#SBATCH --constraint=cascadelake
#SBATCH --partition=commons
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --time=04:00:00

export PATH="/scratch/dg76/bin:$PATH"

THREADS=(1 2 3 4 6 8 12 16 24 32)
FORMATS=("csr" "ell7" "ell8" "tiled")
SIZES=("80 80 80" "100 100 100" "120 120 120" "160 160 160" "200 200 200" "250 250 250")

EXEC_PATH="$(pwd)/test_HPCCG"
export OMP_PLACES=cores
JOB_ID=${SLURM_JOB_ID:-"pid_$$"}

lscpu

for sz in "${SIZES[@]}"; do
    sz_name=$(echo $sz | tr ' ' '_')
    echo "=========================================================="
    echo "MATRIX SIZE: $sz"
    echo "=========================================================="
    
    for fmt in "${FORMATS[@]}"; do
        case $fmt in
            "csr")   fmt_mode=1 ;;
            "ell7")  fmt_mode=2 ;;
            "ell8")  fmt_mode=3 ;;
            "tiled") fmt_mode=4 ;;
            *)       fmt_mode=0 ;;
        esac

        for t in "${THREADS[@]}"; do
            
            export OMP_PROC_BIND=close
            NUMA_CMD="numactl --cpunodebind=0 --membind=0"
            export OMP_NUM_THREADS=$t
            
            PERF_FILE="perf_${fmt}_${sz_name}_t${t}_${JOB_ID}.data"
            
            echo "STARTING: Format=$fmt, Threads=$t, Size=$sz"
            
            perf record -o "$PERF_FILE" -g \
                -e task-clock,cycles,instructions,L1-dcache-load-misses,LLC-load-misses \
                $NUMA_CMD $EXEC_PATH $sz $fmt
            
            if [ -f "$PERF_FILE" ]; then
                python3 parse_perf2.py $fmt_mode "$PERF_FILE"
                rm "$PERF_FILE"
            fi
            
            echo "----------------------------------------------------------"
        done
    done
done