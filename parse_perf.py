import subprocess
import re
import sys
import os

def get_total_duration(data_file):
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
        
        pattern = rf"^\s*([0-9.]+)%\s+([0-9.]+)%\s+([0-9.]+)%\s+([0-9.]+)%\s+([0-9.]+)%\s+([0-9.]+)%.*?({target_prefix}.*)"
        
        found = False
        for line in result.stdout.splitlines():
            match = re.search(pattern, line)
            if match:
                t_clock_pct = float(match.group(1))
                cycles_pct  = float(match.group(2))
                instr_pct   = float(match.group(3))
                l1_refill   = float(match.group(4))
                l2_refill   = float(match.group(5))
                l3_refill   = float(match.group(6))
                symbol      = match.group(7).strip()

                print(f"\nResults for: {symbol}")
                print(f"Source file: {data_file}")
                print("-" * 50)
                print(f"Task-Clock (Time): {t_clock_pct:>7.2f}%")
                print(f"Cycles:            {cycles_pct:>7.2f}%")
                print(f"Instructions:      {instr_pct:>7.2f}%")
                print(f"L1D Refills:       {l1_refill:>7.2f}%")
                print(f"L2D Refills:       {l2_refill:>7.2f}%")
                print(f"L3D Refills:       {l3_refill:>7.2f}%")
                
                if total_seconds:
                    func_time = total_seconds * (t_clock_pct / 100)
                    print(f"Kernel Time:    {func_time:>7.4f} seconds")

                if cycles_pct > 0:
                    print(f"Function IPC:      {instr_pct/cycles_pct:>7.2f}")
                
                if instr_pct > 0:
                    print(f"L1 MPKI:           {(l1_refill/instr_pct)*1000:>7.2f}")
                    print(f"L2 MPKI:           {(l2_refill/instr_pct)*1000:>7.2f}")
                    print(f"L3 MPKI:           {(l3_refill/instr_pct)*1000:>7.2f}")

                found = True
                break
                
        if not found:
            print(f"Function matching '{target_prefix}' not found in the profile.")
        
        return found

    except Exception as e:
        print(f"Error during parsing: {e}")
        return False

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: python3 parse_perf.py [mode] [filename]")
        sys.exit(1)

    mode_map = {
        "1": "csr_spmv", 
        "2": "ellpack7_spmv", 
        "3": "ellpack8_spmv", 
        "4": "ellpack7_tiled_spmv"
    }
    target = mode_map.get(sys.argv[1], "csr_spmv")
    get_perf_stats(target, sys.argv[2])
