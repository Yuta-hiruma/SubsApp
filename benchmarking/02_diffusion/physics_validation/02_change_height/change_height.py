import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import math

# Parameters
z_f = -1.35  # depth of profile [m]
theta_s = 0.44
theta_r = 0.0
n = 1.28
m = 1 - 1 / n
_lambda = 0.1302
kappa = 0.0868
#-----
gamma_s = 2.31
gamma_w = 1
#------
c = 10
phi = 20
ccott=c/math.tan(math.radians(phi))
K0_list = [1, 2]  # m/day
cell = 11
#-------
#dP = 10
dP_list= [1, 5, 10] 
# Depth setup ---------------------------------------------------
z_vals = np.linspace(0, abs(z_f), cell)
z_mids = 0.5 * (z_vals[:-1] + z_vals[1:])
dz = np.diff(z_vals)[0]


#convert parcel size from swap to nsdwl---------------------------

infile_p_ini = "pressure_ini_swp.dat"
infile_wc_ini = "wc_ini_swp.dat"
df_wc_ini = pd.read_csv(infile_wc_ini, delim_whitespace=True, skiprows=1, names=["depth", "wcontent"])
df_pressure_ini = pd.read_csv(infile_p_ini, delim_whitespace=True, skiprows=1, names=["depth", "pressure"])
pressure_depths = df_pressure_ini["depth"].tolist()
wcontent = df_wc_ini ["wcontent"].tolist()
pressures = df_pressure_ini["pressure"].tolist()

theta_ini = []
psi_ini = []
for z in z_vals:
    closest_index = np.argmin(np.abs(abs(np.array(pressure_depths)) - z))
    theta = float(wcontent[closest_index])
    theta_ini.append(theta) 
    psi = pressures[closest_index]
    psi_ini.append(psi)
psi_ini = np.array(psi_ini)
e= np.array(theta_ini)/(1-np.array(theta_ini))


#calculate effective stress and effective stress changes------------
S_sat = (np.array(theta_ini) - theta_r) / (theta_s - theta_r)
sigma_total = 9.81*(gamma_s*(1/(e+1))+S_sat*gamma_w*(1-(1/(e+1))))
sigma_eff = sigma_total - psi_ini


output_dir = "output"
os.makedirs(output_dir, exist_ok=True)


dH = 0
dH_list = []
gwl_list = []

for dP in dP_list:
    psi = psi_ini + dP
    dsigma_eff_dt = np.abs(sigma_eff - (sigma_total - psi))
    dElastic = dsigma_eff_dt * kappa / (ccott + sigma_eff)
    e_results = e - dElastic
    node_heights = (e + 1) / (e_results + 1) * dz
    print(dP)

    for i in range(len(psi)):
        dH = dH + (node_heights[i]-dz)
        if psi[i]<=0:
            #dH = dH + (node_heights[i]-dz)
            gwl =  z_vals[i]
    dH_list.append((dP, dH)) 
    gwl_list.append((dP,gwl))

    # Save profile data
    data = np.column_stack((-z_vals, dElastic, dsigma_eff_dt, sigma_eff, psi,node_heights))
    header = f"Depth[m]\tdElastic\tΔσ_eff\tσ_eff\tpsi\theights\nTotal subsidence above GWT (dH) = {dH:.4f} m"
    filename = os.path.join(output_dir, f"dElastic_profile_dP_{dP}.dat")
    np.savetxt(filename, data, header=header, fmt="%.5f", delimiter="\t")


dH_array = np.array(dH_list)
gwl_array = np.array(gwl_list)

#Plot------------

fig, ax1 = plt.subplots(figsize=(6, 4))

# Primary axis: dH
ax1.set_xlabel('Pressure Increase ΔP [kPa]')
ax1.set_ylabel('Subsidence ΔH [mm]', color='tab:green')
ax1.plot(dH_array[:, 0], dH_array[:, 1]*1000, marker='o', color='tab:green', label='Subsidence dH')
ax1.tick_params(axis='y', labelcolor='tab:green')

# Secondary axis: groundwater level
ax2 = ax1.twinx()
ax2.set_ylabel('Groundwater Level [mm]', color='tab:blue')
ax2.plot(gwl_array[:, 0], -gwl_array[:, 1]*1000, marker='s', linestyle='--', color='tab:blue', label='Groundwater Level')  # use -gwl for below surface
ax2.tick_params(axis='y', labelcolor='tab:blue')

# Title and layout
fig.suptitle('Subsidence and Groundwater Level vs Pressure Change')
fig.tight_layout()
plt.grid(True)
plt.savefig("output/dH_and_gwl_vs_dP.png")