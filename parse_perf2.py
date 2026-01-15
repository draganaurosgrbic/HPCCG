import subprocess
import re
import sys
import os

def get_total_duration(data_file):
    """Extracts the total execution time from perf metadata."""
    try:
        command = ['perf', 'report', '-i', data_file, '--header', '--stdio']
        result = subprocess.run(command, capture_output=True, text=True, check=True)
        
        match = re.search(r"duration\s+:\s+([0-9.]+)\s+ms", result.stdout)
        if match:
            return float(match.group(1)) / 1000.0
    except:
        pass
    return None

def get_perf_stats(target_prefix, data_file):
    if not os.path.exists(data_file):
        print(f"Error: Data file '{data_file}' not found.")
        return False

    command = ['perf', 'report', '-i', data_file, '--stdio', '--group', '--no-children']
    
    try:
        total_seconds = get_total_duration(data_file)
        result = subprocess.run(command, capture_output=True, text=True, check=True)
        pattern = rf"^\s*([0-9.]+)%\s+([0-9.]+)%\s+([0-9.]+)%\s+([0-9.]+)%\s+([0-9.]+)%.*?({target_prefix}.*)"
        
        found = False
        for line in result.stdout.splitlines():
            match = re.search(pattern, line)
            if match:
                t_clock_pct = match.group(1) 
                cycles      = match.group(2) 
                instr       = match.group(3) 
                l1          = match.group(4) 
                llc         = match.group(5) 
                symbol      = match.group(6).strip()

                print(f"\nResults for: {symbol}")
                print(f"Source file: {data_file}")
                print("-" * 50)
                print(f"Task-Clock:    {t_clock_pct:>7}% (Time Share)")
                print(f"Cycles:        {cycles:>7}%")
                print(f"Instructions:  {instr:>7}%")
                print(f"L1 Misses:     {l1:>7}%")
                print(f"LLC Misses:    {llc:>7}%")
                
                if total_seconds:
                    kernel_time = total_seconds * (float(t_clock_pct) / 100)
                    print(f"Kernel Time:   {kernel_time:>7.4f} seconds")
                    print(f"Total Runtime: {total_seconds:>7.4f} seconds")

                try:
                    c_val = float(cycles)
                    t_val = float(t_clock_pct)
                    if t_val > 0:
                        print(f"Cycle/Clock Intensity: {c_val/t_val:>7.2f}")
                    if c_val > 0:
                        print(f"Relative IPC:          {float(instr)/c_val:>7.2f}")
                except:
                    pass
                found = True
                break
                
        if not found:
            print(f"Function matching '{target_prefix}' not found.")
            if total_seconds:
                print(f"Total Runtime: {total_seconds:>7.4f} seconds (Kernel not found)")
        
        return found

    except Exception as e:
        print(f"Error: {e}")
        return False

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: python3 parse_perf2.py [mode] [filename]")
        sys.exit(1)

    mode_map = {"1": "csr_spmv", "2": "ellpack7_spmv", "3": "ellpack8_spmv", "4": "ellpack7_tiled_spmv"}
    target = mode_map.get(sys.argv[1], "csr_spmv")
    get_perf_stats(target, sys.argv[2])
