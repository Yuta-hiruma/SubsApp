import sys
import os
current_dir = os.getcwd()
parent_dir = os.path.abspath("..")
sys.path.insert(0, parent_dir)
from Input_soil_material_parameters import *  # Assumes 'depths' is defined here
from input import *                           # Assumes 'swap_r.swp' is defined here
sys.path.pop(0)
os.chdir(current_dir)


heights = [(-depths[0]) * 100]  # First sublayer thickness from 0 to first depth
for i in range(len(depths) - 1):
    heights.append(round((depths[i] - depths[i+1]) * 100, 3))


#Simulation period and site name--------------------------------------------------------
with open('swap.swp', 'r') as f:
    lines = f.readlines()
output_lines = []

for line in lines:
    stripped = line.strip()
    
    #TSTART replacement---------------------------------------
    if stripped.startswith("TSTART ="):
        parts = line.split('=')
        if len(parts) > 1:
            comment = parts[1].split('!')[1] if '!' in parts[1] else ''
            new_line = f"{parts[0]}= {tstart}    !{comment}\n"
            output_lines.append(new_line)
        else:
            output_lines.append(line)
            
    #TEND replacement---------------------------------------
    elif stripped.startswith("TEND ="):
        parts = line.split('=')
        if len(parts) > 1:
            comment = parts[1].split('!')[1] if '!' in parts[1] else ''
            new_line = f"{parts[0]}= {tend}    !{comment}\n"
            output_lines.append(new_line)
        else:
            output_lines.append(line)

    #Site name---------------------------------------
    elif stripped.startswith("METFIL ="):
        parts = line.split('=')
        if len(parts) > 1:
            comment = parts[1].split('!')[1] if '!' in parts[1] else ''
            new_line = f"{parts[0]}= '{site}'   !{comment}\n"
            output_lines.append(new_line)
        else:
            output_lines.append(line)
            
    #-----------------------------------------------
    elif stripped.startswith("RAINFIL ="):
        parts = line.split('=')
        if len(parts) > 1:
            comment = parts[1].split('!')[1] if '!' in parts[1] else ''
            new_line = f"{parts[0]}= '{site}'    !{comment}\n"
            output_lines.append(new_line)

    #-----------------------------------------------
    elif stripped.startswith("BBCFIL ="):
        parts = line.split('=')
        if len(parts) > 1:
            comment = parts[1].split('!')[1] if '!' in parts[1] else ''
            new_line = f"{parts[0]}= '{site}'    !{comment}\n"
            output_lines.append(new_line)

    #-----------------------------------------------
    elif stripped.startswith("GWLI   ="):
        parts = line.split('=')
        if len(parts) > 1:
            comment = parts[1].split('!')[1] if '!' in parts[1] else ''
            new_line = f"{parts[0]}= {gwli}    !{comment}\n"
            output_lines.append(new_line)

    else:
        output_lines.append(line)
with open('swap.swp', 'w') as f:
    f.writelines(output_lines)



#Part 4: Vertical discretization of soil profile------------------------------
with open('swap.swp', 'r') as f:
    lines = f.readlines()

output_lines = []
in_table = False
j = 0
for line in lines:
    stripped = line.strip()

    if stripped.startswith("ISUBLAY") and "HSUBLAY" in stripped:
        in_table = True
        output_lines.append(line)
        continue

    if in_table and stripped.startswith("* End of table"):
        in_table = False
        output_lines.append(line)
        continue

    if in_table and not stripped.startswith("*") and stripped:
        soil_profile = line.split()
        if  j < len(heights):
            hsub = heights[j]
            hcomp = 5.0
            ncomp = int(round(hsub / hcomp))
            soil_profile[0] = j+1
            soil_profile[1] = j+1
            soil_profile[2] = f"{hsub:.1f}"   # HSUBLAY
            soil_profile[3] = f"{hcomp:.1f}"  # HCOMP
            soil_profile[4] = str(ncomp)      # NCOMP

            j += 1
            new_line = "       "+"{:<10}{:<5}{:<8}{:<8}{:<}".format(soil_profile[0], soil_profile[1], soil_profile[2], soil_profile[3], soil_profile[4]) + "\n"
            output_lines.append(new_line)
        else:
            output_lines.append(line)
    else:
        output_lines.append(line)

with open('swap.swp', 'w') as f:
    f.writelines(output_lines)


#Part 5: Soil hydraulic functions--------------------------------------------
with open('swap.swp', 'r') as f:
    lines = f.readlines()

output_lines = []
in_table = False
j = 0
for line in lines:
    stripped = line.strip()
    if stripped.startswith("ORES") and "OSAT" in stripped:
        in_table = True
        output_lines.append(line)
        continue

    if in_table and stripped.startswith("* End of table"):
        in_table = False
        output_lines.append(line)
        continue

    if in_table and not stripped.startswith("*") and stripped:
        hyd= line.split()
        if  j < len(heights):
            soil_type = types[j]
            hyd[0] = Soil_parameters[soil_type]["Water_retention_curve"]["vG_WCR"] #ORES
            hyd[1] = Soil_parameters[soil_type]["Water_retention_curve"]["vG_WCS"] #OSAT
            hyd[2] = Soil_parameters[soil_type]["Water_retention_curve"]["vG_alpha"] #ALFA
            hyd[3] = Soil_parameters[soil_type]["Water_retention_curve"]["vG_n"] #NPAR
            hyd[4] = Soil_parameters[soil_type]["Water_retention_curve"]["Ksatfit"] #KSATFIT
            hyd[5] = Soil_parameters[soil_type]["Water_retention_curve"]["Lexp"] #LEXP
            hyd[6] = Soil_parameters[soil_type]["Water_retention_curve"]["H_enpr"] #H_ENPR
            hyd[7] = Soil_parameters[soil_type]["Water_retention_curve"]["Ksattexm"] #KSAREXM
            hyd[8] = Soil_parameters[soil_type]["Water_retention_curve"]["Bdens"] #BDENS
            j += 1
            new_line = "  "+"{:<6}{:<7}{:<8}{:<8}{:<9}{:<9}{:<10}{:<8}{:<}".format(hyd[0], hyd[1], hyd[2], hyd[3], hyd[4], hyd[5], hyd[6], hyd[7], hyd[8]) + "\n"
            output_lines.append(new_line)
        else:
            output_lines.append(line)
    else:
        output_lines.append(line)

with open('swap.swp', 'w') as f:
    f.writelines(output_lines)
#----------------------------------------------------------------------------
