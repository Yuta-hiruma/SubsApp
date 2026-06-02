import pandas as pd
import re
from datetime import datetime

file_path = "../result.vap"
swp_path = "../swap.swp"
out_file = "../../nsdwl/input/wc_ini/wc_ini_swp.dat"
out_file_ini = "../../../src/input/wc_ini/wc_ini_ave_swp.dat"   #convert output of swap to input for nsdwl

hlayer = []
dlayer = []
skip = 0
read_data = False

with open(swp_path, 'r') as swp:
    for line in swp:
        line = line.strip() 
        if "TSTART " in line:
            TSTART = line.split('=')[1].split('!')[0].strip().strip("'")
        elif "TEND " in line:
            TEND = line.split('=')[1].split('!')[0].strip().strip("'")
        if re.search(r"\bISUBLAY\b", line):
            read_data = True 
            continue  

        if read_data:
            if line.startswith("* End of table"):
                break  
            parts = re.split(r'\s+', line)
            if len(parts) >= 3:  
                try:
                    hsub_value = float(parts[2])
                    hlayer.append(hsub_value)
                except ValueError:
                    continue

TSTART_date = datetime.strptime(TSTART, "%Y-%m-%d")
TEND_date = datetime.strptime(TEND, "%Y-%m-%d")

for i in range(len(hlayer)):
    if i == 0:
        dlayer.append(- hlayer[i])
    else:
        dlayer.append(dlayer[i-1] - hlayer[i])
date_plt = []
phead_plt = []
ave_wc = []
n_layer = 0
n = 0
sum_wc = 0.0
with open(file_path, "r", encoding="utf-8", errors="replace") as infile, open(out_file, 'w') as outfile, open(out_file_ini, 'w') as outfile_ini:
    lines = infile.readlines()[12:]  
    outfile.write(f"{'Depth(m)':<11} {'water content'}\n") 
    outfile_ini.write(f"{'Layer':<11} {'water content'}\n") 

    for row in lines:
        values = row.split(',')  
        if len(values) < 4:
            continue 

        date = values[0].strip()
        date_obj = datetime.strptime(date, "%Y-%m-%d")
        depth = float(values[1]) /100
        wcontent = values[2]
        phead = values[3]
        if not date.strip():
            continue  
        if date_obj == TSTART_date:
            outfile.write(f"{depth:<11} {wcontent}\n") 
            if depth >= dlayer[n_layer]:
                if any(v.strip() == "" for v in values):
                    ave_wc.append(sum_wc/n)
                    n_layer = n_layer + 1
                    outfile_ini.write(f"{n_layer:<11} {ave_wc[n_layer-1]}\n")
                else:
                    n = n+1
                    sum_wc += float(wcontent)
            else:
                ave_wc.append(sum_wc/n)
                n_layer = n_layer + 1
                outfile_ini.write(f"{n_layer:<11} {ave_wc[n_layer-1]}\n") 
                n = 1
                sum_wc = float(wcontent)







