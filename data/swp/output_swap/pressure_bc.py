import pandas as pd
import matplotlib.pyplot as plt
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
out_file_tp = "fig/pressure_head.dat"
out_file_tp_swp = "../../nsdwl/input/top_pressure/top_pressure_swp.dat"   #convert output of swap to input for nsdwl

#out_file_bp = "fig/pressure_head.dat"
#out_file_bp_swp = "../../nsdwl/input/bottom_pressure/bottom_pressure_swp.dat" 


rho = 1000 #kg/m3
g = 9.81 #m/s2

date_plt = []
phead_plt = []
date_count = 0

#top boundary pressure------------------------------------------
with open(file_path, "r", encoding="utf-8", errors="replace") as infile, open(out_file_tp, 'w') as outfile, open(out_file_tp_swp, 'w') as outfile_tp:
    lines = infile.readlines()[12:]  
    outfile_tp.write(f"{'date':<11} {'pressrue at top(kPa)'}\n") 
    for row in lines:
        values = row.split(',')  

        if len(values) < 4:
            continue 

        date = datetime.strptime(values[0].strip(), "%Y-%m-%d")
        depth = float(values[1])
        if date < tstart:
            continue


        if values[3].strip():
            try:
                phead = float(values[3].strip()) / 100.0 * rho * g /1000#convert to kPa
            except ValueError:
                continue
        else:
            continue


        #if not date.strip():
        #    continue  
        if depth == -2.5:
            outfile.write(f"{date:<11} {phead}\n") 
            outfile_tp.write(f"{date_count:<11} {round(float(phead),2)}\n")
            date_plt.append(date)
            phead_plt.append(float(phead))
            date_count = date_count + 1


# Plot Pressure Head
plt.figure(figsize=(15, 6))
plt.plot(date_plt, phead_plt, linestyle="-")
plt.xlabel("Date")
plt.ylabel("Pressure (kPa)")
plt.title("Pressure at Depth = -2.5 cm")
plt.xticks(rotation=45)
plt.xticks(date_plt[::60])
plt.subplots_adjust(left=0.2, bottom=0.2)
plt.legend(["Pressure Head"])
#plt.grid(True)
plt.savefig("fig/phead_top.png", dpi=300)

date_plt_bt = []
phead_plt_bt = []
date_count_bt = 0

#bottom boundary pressure------------------------------------------
#with open(file_path, "r", encoding="utf-8", errors="replace") as infile, open(out_file_bp, 'w') as outfile, open(out_file_bp_swp, 'w') as outfile_bp:
#    lines = infile.readlines()[12:]  
#    outfile_bp.write(f"{'date':<11} {'pressrue at bottom(kPa)'}\n") 
#    for row in lines:
#        values = row.split(',')  
#
#        if len(values) < 4:
#            continue 
#
#        date = values[0]
#        depth = float(values[1])
#        if values[3].strip():
#            try:
#                phead = float(values[3].strip()) / 100.0 * rho * g /1000#convert to kPa
#            except ValueError:
#                continue
#        else:
#            continue
#
#
#        if not date.strip():
#            continue  
#        if depth == -582.5:
#            outfile.write(f"{date:<11} {phead}\n") 
#            outfile_bp.write(f"{date_count_bt:<11} {round(float(phead),2)}\n")
#            date_plt_bt.append(date)
#            phead_plt_bt.append(float(phead))
#            date_count_bt = date_count_bt + 1
#
#
## Plot Pressure Head
#plt.figure(figsize=(15, 6))
#plt.plot(date_plt_bt, phead_plt_bt, linestyle="-")
#plt.xlabel("Date")
#plt.ylabel("Pressure Head (cm)")
#plt.title("Pressure Head at Depth = -582.5 cm")
#plt.xticks(rotation=45)
#plt.xticks(date_plt_bt[::60])
#plt.subplots_adjust(left=0.2, bottom=0.2)
#plt.legend(["Pressure Head"])
##plt.grid(True)
#plt.savefig("fig/phead_bottom.png", dpi=300)
