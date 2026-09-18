import pandas as pd
import re
from datetime import datetime

#------
import sys
import os
current_dir = os.getcwd()
parent_dir = os.path.abspath("../..")
sys.path.insert(0, parent_dir)
from Input_soil_material_parameters import *
from input import *
grandparent_dir = os.path.abspath(os.path.join("../../.."))
sys.path.insert(0, grandparent_dir)
from general_input import *  
sys.path.pop(0)
os.chdir(current_dir)
#------
tstart = datetime.strptime(tstart, "%Y-%m-%d")

file_path = "../result.vap"
swp_path = "../swap.swp"
out_file = "../../nsdwl/input/pressure_ini/pressure_ini_swp.dat"  #convert output of swap to input for nsdwl

rho = 1000 #kg/m3
g = 9.81 #m/s2

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
with open(file_path, "r", encoding="utf-8", errors="replace") as infile, open(out_file, 'w') as outfile:
    lines = infile.readlines()[12:]  
    outfile.write(f"{'Depth(m)':<11} {'pressure'}\n") 

    for row in lines:
        values = row.split(',')  
        if len(values) < 4:
            continue 

        date = values[0].strip()
        date_obj = datetime.strptime(date, "%Y-%m-%d")
        depth = float(values[1]) /100

        if values[3].strip():
            try:
                pressure = float(values[3].strip()) / 100.0 * rho * g /1000#convert to kPa
            except ValueError:
                continue
        else:
            continue


        if not date.strip():
            continue  
        if date_obj == TSTART_date:
            outfile.write(f"{depth:<11} {round(float(pressure),4)}\n") 








