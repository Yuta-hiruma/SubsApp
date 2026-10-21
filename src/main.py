 # -*- coding: utf-8 -*-
"""
Created on Fri Jul 23 14:27:28 2021  
@author: tomdegast


Updated on Fri Jun 11 2026 
@editor: Yuta_hiruma

"""

#_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
#Preparation ----------------
from scipy import special
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as ticker 
import matplotlib.colors as mcolors
from matplotlib.colors import ListedColormap, BoundaryNorm
import numpy as np
import scipy as sci
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

from pathlib import Path
import math as math
import pandas as pd
import time
from datetime import datetime, timedelta
from dask import delayed
import dask.bag as db
#Numba ----------------------
from numba import int32, float32
from numba import jit
#Local files ----------------
from build_domain import *
from NSDWLs_errors import *
from NSDWLs_equations import *
#Reading input data ----------
import sys
import os
current_dir = os.getcwd()
parent_dir = os.path.abspath("../data")
sys.path.insert(0, parent_dir)
from Input_soil_material_parameters import *
from input import *
sys.path.pop(0)
os.chdir(current_dir)
#Reading boundary condition files-----------------------------
infile_tp = "input/top_pressure/top_pressure_swp.dat"
infile_bp = "input/bottom_pressure/bottom_pressure_swp.dat"
infile_wc_ini = "input/wc_ini/wc_ini_swp.dat"
infile_p_ini = "input/pressure_ini/pressure_ini_swp.dat"
df_tp = pd.read_csv(infile_tp, sep=r"\s+", skiprows=1, names=["date", "pressure"])
df_bp = pd.read_csv(infile_bp, sep=r"\s+", skiprows=1, names=["date", "pressure"])
df_wc_ini = pd.read_csv(infile_wc_ini, sep=r"\s+", skiprows=1, names=["depth", "wcontent"])
df_pressure_ini = pd.read_csv(infile_p_ini, sep=r"\s+", skiprows=1, names=["depth", "pressure"])
#_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

#Initilization of the soil domain-------------------------------------
Soil_dimension = Soil_domain()
Soil_dimension.Depth_top = 0.00                     # m size
Soil_dimension.Depth_bottom = depths[-1]            # m size

"Default soil column discretization is 0.1 m. Adjust this value according to the required model resolution."
Soil_dimension.precision = 0.05 # m node size
#------------------------------
Soil_dimension.Depth_o_soils = depths                # bottom of soil layer      
Soil_dimension.soils = types                                      # soil types used
Soil_dimension.Number_o_soils = len(Soil_dimension.Depth_o_soils)                                                 # number of soils 
Soil_dimension.domainsize = Soil_dimension.Depth_top-Soil_dimension.Depth_bottom

Soil_states = Soil_states()
Soil_stresses = Soil_stresses()

Soil_calc = Soil_calc()
swp = swp()

#Simulaiton time ------------------------------------------------
"Aduest max_time_step and relative tolerance that are used in differential equation." 
max_time_step = 1.0  #day
rtol = 1e-07#relative tolerance 
#atol =1e-4 #absolute tolerance
#--------------------------------
TSTART = tstart
TEND = tend

start_date = datetime.strptime(TSTART, "%Y-%m-%d")
end_date = datetime.strptime(TEND, "%Y-%m-%d")
end_time = (end_date - start_date).days
simssteps = (end_date - start_date).days
date_list = [start_date + timedelta(days=i) for i in range((end_date - start_date).days + 1)]
simulation_time =  np.linspace (0,end_time,simssteps+1)

#Set boundary conditions ------------------------------------------
#Top boundary conditions ------------------------------
Bctop = Bc_top()
#Initial conditoin from SWAP----------
Bctop.time = df_tp["date"].tolist()
Bctop.pressure = df_tp["pressure"].tolist()   
#-------------------------------------
Bctop.trans = np.ones(len(Bctop.pressure))
Bctop.load_increments = [0]
Bctop.load = np.cumsum(Bctop.load_increments)
Bctop.time_incr = [0]
#Bottom boundary conditions ----------------------------
Bctop.devaporation = np.zeros(len(Bctop.pressure))
Bcbottom = Bc_bottom()
#Bcbottom.time = [0]
#Bcbottom.pressure = [49]
Bcbottom.time = df_bp["date"].tolist()
Bcbottom.pressure = df_bp["pressure"].tolist()     # kPa

bc = 0 # if Boudary condition is given as time series bc = 1
dpsi_top_fun=[]
dpsi_bottom_fun=[]
if len(Bcbottom.time) != 1:
    bc = 1

if bc ==1:
    #Boundary-condition interpolation
    psi_top_fun = CubicSpline(Bctop.time, Bctop.pressure)
    psi_bottom_fun = CubicSpline(Bcbottom.time, Bcbottom.pressure)
    
    #Time derivatives of boundary condition
    dpsi_top_fun = psi_top_fun.derivative()
    dpsi_bottom_fun = psi_bottom_fun.derivative()


#Check errors -----------------------------------------------------
checkerror(Soil_dimension,Soil_parameters)

#Set initial conditions -------------------------------------------
#Initial conditoin from SWAP----------
swp.depth = df_wc_ini ["depth"].tolist()
swp.wcontent = df_wc_ini ["wcontent"].tolist()
swp.pressure = df_pressure_ini ["pressure"].tolist() #kPa
#Initail Soil states and domain--------
(Soil_states,Soil_states.theta_ini,Soil_stresses.psi_ini) = set_ini(Soil_dimension,Soil_parameters,Soil_states,Bctop,swp)
Soil_states.e_ini = Soil_states.e
Soil_states.n_ini = Soil_states.n
Soil_states.theta = Soil_states.theta_ini

Soil_dimension.internode_vector = modeldomain (Soil_dimension,Soil_parameters)

(material_vector,gamma_vector,Permeability_curve_Perm_vector,Water_retention_curve_vector, Compression_vector) = soilbuilder(Soil_dimension,Soil_parameters,Soil_dimension.node_vector,Soil_dimension.internode_vector)

S_e_ini = Soil_states.theta/Soil_states.n_ini
Soil_dimension.node_initial=get_node_distance(Soil_dimension)
(Soil_stresses) = calc_stress_ini_h (Soil_dimension,gamma_vector,Soil_states,Water_retention_curve_vector,Soil_stresses)
Soil_stresses.sigma_effective_node_ini = np.array(Soil_stresses.sigma_effective_node)
Soil_stresses.sigma_effective_p = Compression_vector.OCR*Soil_stresses.sigma_effective_node_ini
sigma_eff_p=Soil_stresses.sigma_effective_p

#Define a wrapper function to store previous values------
class SigmaTracker:
    def __init__(self):
        self.prev_sigma_eff = None 
        self.prev_t = None     
    
    def compute_dsigma_eff_dt(self, sigma_eff, t):
        if self.prev_sigma_eff is None:  
            dsigma_eff_dt = 0  
        else:
            dt = t - self.prev_t 
            if dt == 0:
                dsigma_eff_dt = 0 
            else:
                dsigma_eff_dt = (sigma_eff - self.prev_sigma_eff) / dt 
        self.prev_sigma_eff = sigma_eff
        self.prev_t = t
        return dsigma_eff_dt
#Initialize tracker ------------------------------
sigma_tracker = SigmaTracker()


#Defining the equation for PDE ------------------------
def equations(simulation_time, in_in, sigma_tracker):
    """
        Computes the time derivatives of pore pressure (psi), water content (theta),
        and void ratio (e) for all soil column elements.
        
        Parameters
        ----------
        simulation_time : float
            Current simulation time provided by the ODE solver. The simulation time
            vector is  defined as np.linspace(0, end_time, sim_steps + 1)
        
        in_in : ndarray
            State vector containing the values of pore pressure (psi),
            water content (theta), and void ratio (e) for each soil element.
        
        sigma_tracker : object
            Tracks changes in effective stress (sigma_eff). Effective stress is not
            solved directly as a state variable in this system of equations. Instead,
            sigma_eff is updated after each time step based on the updated values of
            effective saturation (S_e), void ratio (e), and pore pressure (psi).
            The rate of change of effective stress is then obtained from the
            difference between the previous and updated sigma_eff values.
        
        Returns
        -------
        ndarray
            Time derivatives of the state variables (psi, theta, and e) for all
            soil elements.
    """    

    global gamma_vector
    global Compression_vector
    global Water_retention_curve_vector
    global Permeability_cruve_Perm_vector
    global Soil_dimension
    global Bctop
    global Bcbottom
    print ("time", "%.2f" % simulation_time)
    #setup_soil_states  ------------ 
    Water_retention_curve_node = Water_retention_curve_vector
    nr_elements = int(np.size(in_in)/6)  
    Pressure,theta,e = untangle_parameters(in_in.copy(),3)
    Soil_states.theta = theta
    Soil_states.e = e
    n = e/(1+e)
    S_e = np.clip(update_S_e(Water_retention_curve_node, theta), 0.0, 1.0)
    saturated = (S_e >= 1).astype(int)
    unsaturated = (S_e < 1).astype(int)  
    #_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/  
    #Top pressure boundary condition------------
    t_interval = len(Bctop.time)-1
    for i in range(len(Bctop.time)):
        if simulation_time < Bctop.time[i]:
            t_interval = t_interval-1   
    Bctop_in = Bc_top()
    Bctop_in.pressure =  Bctop.pressure[t_interval]
    if t_interval == len(Bctop.time)-1 :
        t_next = Bctop.time[t_interval]
        P_next = Bctop.pressure[t_interval]
        trans_next = Bctop.trans[t_interval]
    else: 
        t_next = Bctop.time[t_interval+1]
        P_next = Bctop.pressure[t_interval+1]
        trans_next = Bctop.trans[t_interval+1]
        P_average = 0.5*(Bctop.pressure[t_interval] + P_next)
        P_range = P_next - Bctop.pressure[t_interval]
        Bctop_in.pressure = P_average+0.5*P_range*special.erf((simulation_time-t_next+trans_next)*(3/trans_next))
    Bctop_in.devaporation = Bctop.devaporation[t_interval]
    Bctop_in.trans = Bctop.trans[t_interval]
    #Bottom pressure boundary condition------------
    t_interval = len(Bcbottom.time)-1
    for i in range(len(Bcbottom.time)):
        if simulation_time < Bcbottom.time[i]:
            t_interval = t_interval-1  
    Bcbottom_in = Bc_bottom()
    Bcbottom_in.pressure = Bcbottom.pressure[t_interval] 

 
    #Update soil_dimension------------
    Soil_dimension.node_vector = get_node_coordinates(Soil_dimension)
    Soil_dimension.node_distance = get_node_distance(Soil_dimension)  
    Soil_dimension.node_distance = (1-((Soil_states.e_ini-e)/(1+Soil_states.e_ini)))*Soil_dimension.node_initial
 
    #Update_soil_stresses------------
    sigma_total = np.zeros(len(e))
    sigma_total = update_sigma_total (Soil_dimension.node_distance,gamma_vector,Soil_states,Water_retention_curve_vector,Bctop_in)
    sigma_eff = np.abs(sigma_total - Pressure*S_e)
    #compute_dsigma_eff_dt------------------------------------------------
    dsigma_eff_dt = sigma_tracker.compute_dsigma_eff_dt(sigma_eff, simulation_time) 
    #---------------------------------------------------------------------
    psi_m = np.where(Pressure <= 0, Pressure, 0)

    #Calculate permeability----------
    Soil_calc.K                     = update_K_Kmin (theta,e,S_e,Permeability_curve_Perm_vector,Water_retention_curve_node)
    K= Soil_calc.K

    #set differentials ------------------------
    #First set: pressure ----------
    #d2psi_dz2 = d2f_dz2_non_uniform(Pressure,Bctop_in.pressure,Bcbottom_in.pressure,Soil_dimension) #Pressure*0+np.sin(simulation_time)
    #dpsi_dt = K*d2psi_dz2
    #Diffusion equation
    dpsi_dt = diffusion(Pressure,Bctop_in.pressure,Bcbottom_in.pressure,dpsi_top_fun,dpsi_bottom_fun,bc,Soil_dimension,K,simulation_time)

    dS_e_dpsi,dtheta_dpsi  = derrivative_van_genuchten_S_e(Water_retention_curve_node, theta,e, psi_m)
    # second set: theta ----------
    dtheta_dt = dpsi_dt*dtheta_dpsi  
    # third set: void ratio -------
    for i in range(len(sigma_eff_p)):
        if sigma_eff[i] > sigma_eff_p[i]:
            sigma_eff_p[i] = sigma_eff[i]  
    de_dt,dn_dt = deSSc_dt(Compression_vector, Soil_states, S_e, sigma_eff,sigma_eff_p, dsigma_eff_dt,dpsi_dt,Water_retention_curve_vector)


    compression = (de_dt < 0).astype(int)    
    extension = (de_dt > 0).astype(int)  
    
    positive_pressure = (dpsi_dt > 0).astype(int)
    negative_pressure = (dpsi_dt < 0).astype(int)

    dpsi_dt = update_dpsi_dt(dpsi_dt, K, compression, saturated, dn_dt)

    dtheta_dt, dSat_dt = update_dtheta_and_dSat(dtheta_dt,positive_pressure, extension, saturated, dn_dt,n, theta, unsaturated)
    dtheta_dt =  gather(*(saturated_theta(S_e,dtheta_dt,dn_dt,i) for i in range(len(S_e))))#.compute()

    out_out = np.concatenate([dpsi_dt,dtheta_dt,de_dt],axis=0)
    return [out_out]



#Set input ---------
# First set is pressure #
input_in = Soil_stresses.psi_ini
# Second set is water content #
input_in = np.append(input_in,Soil_states.theta)
# Third set is void ratio #
input_in = np.append(input_in,Soil_states.e_ini) 


# Solveing differential equation ------
start = time.perf_counter() # measure the performance
sol = solve_ivp(equations, [simulation_time[0], simulation_time[-1]],input_in,method='RK23',t_eval=simulation_time,max_step = max_time_step, rtol = rtol ,args=(sigma_tracker,))#, atol = atol)
#sol = solve_ivp(equations, [simulation_time[0], simulation_time[-1]],input_in,method='Radau',t_eval=simulation_time,max_step = max_time_step, rtol = rtol ,args=(sigma_tracker,))#, atol = atol)
end = time.perf_counter()# measure the performance

print("Elapsed  = {}s".format((end - start)))
print("Success:", sol.success)
print("Message:", sol.message)
print("Final time reached:", sol.t[-1])

#Post process-------------------------------------------
nr_elements = int(np.size(input_in)/3)
# First set is pressure -----
Pressure_results = sol.y[0:nr_elements]
# Second set is water content-----
theta_results = sol.y[nr_elements:2*nr_elements]
# Third set is void ratio-----
e_results = sol.y[2*nr_elements:3*nr_elements]
#-------------------------------------------------------
n_results = e_results/(1+e_results)
Se_results = theta_results/n_results
MR = theta_results*(1+e_results) #Moisture ratio
#-------------------------------------------------------

#Convert void ratio to surface displacement -----
colum_height=np.zeros(len(e_results[0,:]))
node_height=0*(e_results)
for i in range (len(e_results[0,:])):
    colum_height[i] = np.sum((1-((Soil_states.e_ini-e_results[:,i])/(1+Soil_states.e_ini)))*Soil_dimension.node_initial)
    node_height[:,i] = ((1-((Soil_states.e_ini-e_results[:,i])/(1+Soil_states.e_ini)))*Soil_dimension.node_initial)
colum_node_height = node_height*0
colum_node_height[-1,:] = node_height[-1,:]

for i in range(len(e_results[:,0])-2,-1,-1):
    colum_node_height[i,:] = node_height[i,:]+colum_node_height[i+1,:]

#ground water level-------------------------------------------
gwl = []

for i in range(Pressure_results.shape[1]):  # Loop over columns (time)
    for j in range(Pressure_results.shape[0] - 1):  # Loop over rows (depth)
        if Pressure_results[j, i] > 0:
            gwl.append(colum_node_height[j, i] + Soil_dimension.Depth_bottom)
            break 

#-------------------------------------------------------------
#Output 
#-------------------------------------------------------------
#dat file-----------------------------------------------------
with open('../results/displacement.dat', 'w') as dis:
    dis.write(f"{'':<2}{'time[days]':<15}  {'displacement [mm]':<15}\n")
    for i in range(len(date_list)):
        dis.write(f"{'':<2}{date_list[i].strftime('%Y-%m-%d'):<15}  {(colum_height[i]+Soil_dimension.Depth_bottom)*1000:<15}\n")

with open("../results/node_heights.dat", "w") as f:
    f.write(f"{'time[day]':<15}")
    for j in range(node_height.shape[0]):
        f.write(f"{f'Node_{j}[mm]':<15}")
    f.write("\n")
    for i in range(len(simulation_time)):
        f.write(f"{simulation_time[i]:<15.6f}")
        for j in range(node_height.shape[0]):
            f.write(f"{node_height[j, i] * 1000:<15.6f}")
        f.write("\n")

with open('../results/Pressure_top.dat', 'w') as dis:
    dis.write(f"{'':<2}{'time [days]':<15} {'Pressure top':<15} \n")
    for i in range(len(date_list)):
        dis.write(f"{'':<2}{date_list[i].strftime('%Y-%m-%d'):<15} {Pressure_results[0,i]:<15} \n")

with open('../results/Pressure_bottom.dat', 'w') as dis:
    dis.write(f"{'':<2}{'time [days]':<15} {'Pressure top':<15} \n")
    for i in range(len(date_list)):
        dis.write(f"{'':<2}{date_list[i].strftime('%Y-%m-%d'):<15} {Pressure_results[-1,i]:<15} \n")

with open('../results/n.dat', 'w') as dis:
    dis.write(f"{'':<2}{'time [days]':<15}  {'depth':<15} {'n':<15}\n")
    for i in range(len(date_list)):
        dis.write(f"{'':<2}{date_list[i].strftime('%Y-%m-%d'):<15} {colum_node_height[1,i]+Soil_dimension.Depth_bottom:<15} {e_results[1,i]/(1+e_results[1,i]):<15}\n")


#Plot displacement---------------------------------- (This should be moved to /gmt)
plt.figure()
#plt.plot(date_list,(colum_height+Soil_dimension.Depth_bottom)*1000)
plt.plot(simulation_time,(colum_height+Soil_dimension.Depth_bottom)*1000)
plt.xlabel('time [days]')
plt.ylabel('displacement [mm]')
plt.title('displacement')
#plt.xticks(rotation=45)
plt.subplots_adjust(left=0.2, bottom=0.2)
plt.savefig('../results/figure/displacement.png')


#Plot preparation------------------------------------
time = np.asarray(simulation_time)
vG_WCR = Water_retention_curve_vector.vG_WCR[:, np.newaxis]  
vG_WCS = Water_retention_curve_vector.vG_WCS[:, np.newaxis]
depth = np.round(np.asarray(Soil_dimension.node_vector),8)
effective_saturation = np.round(((theta_results - vG_WCR) / (vG_WCS - vG_WCR)),8)
pressure = np.round(np.asarray(Pressure_results),8)
theta = np.round(np.asarray(theta_results),8)

##Hydraulic conductivity_/_/ _/_/ _/_/ _/_/ _/_/ _/_/ _/_/ _/_/ 
e = np.round(np.asarray(e_results),8)
K_sat = Permeability_curve_Perm_vector.K_tot[:, np.newaxis]
K_min = Permeability_curve_Perm_vector.K_min[:, np.newaxis]
tau = Permeability_curve_Perm_vector.tau[:, np.newaxis] 
vG_m = Water_retention_curve_vector.vG_m[:, np.newaxis]
K_out = K_sat*(effective_saturation**tau)*(1-(1-effective_saturation**(1/vG_m)**vG_m))**2
K_out = np.maximum(K_out, K_min)
K_out = np.where(effective_saturation < 1,K_out,K_sat)
#_/_/ _/_/ _/_/ _/_/ _/_/ _/_/ _/_/ _/_/ _/_/ _/_/ _/_/ _/_/ 



time_grid, depth_grid = np.meshgrid(time, depth)
# Pressure------------------------------------
data = np.column_stack([time_grid.ravel(),depth_grid.ravel(),pressure.ravel(),])
np.savetxt("../results/Pressure.dat",data,delimiter="\t",header="time_days\tdepth_m\tpressure_kPa",fmt="%.8e",)
# Se_vG--------------------------------------
data = np.column_stack([time_grid.ravel(),depth_grid.ravel(),np.asarray(effective_saturation).ravel()])
np.savetxt(f"../results/Se_vG.dat",data,delimiter="\t",header="time\tdepth\teffective_saturation",fmt="%.4e")
# Theta--------------------------------------
data = np.column_stack([time_grid.ravel(),depth_grid.ravel(),np.asarray(theta).ravel()])
np.savetxt(f"../results/theta.dat",data,delimiter="\t",header="time\tdepth\ttheta",fmt="%.4e")
# Hydraulic conductivity----------------------
data = np.column_stack([time_grid.ravel(),depth_grid.ravel(),np.asarray(K_out).ravel()])
np.savetxt(f"../results/K.dat",data,delimiter="\t",header="time\tdepth\tK",fmt="%.4e")

#--------------------------------
#Benchmarking
#--------------------------------
i = 2
# --------------------------------------------------
# Stress - strain
sigma_total_node = 9.81 * (gamma_vector.s[:, None] * (1 / (e_results + 1))+ Se_results * gamma_vector.w[:, None] * (1 - (1 / (e_results + 1))))
sigma_total_in = np.zeros((i+2,len(sigma_total_node[0])))
sigma_total = np.zeros((i+1,len(sigma_total_node[0])))
for j in range(i+1):
    sigma_total_in[j+1] = sigma_total_node[j]*Soil_dimension.node_distance[j]+sigma_total_in[j]
    sigma_total[j] = 0.5*sigma_total_in[j]+0.5*sigma_total_in[j+1]

sigma_eff = np.abs(sigma_total[i] - Pressure_results[i]*Se_results[i])
strain = (Soil_states.e_ini[i]-e_results[i])/Soil_states.e_ini[i]

# --------------------------------------------------
# Soil shrinkage characteristic curve - Se from theta
Se_theta = ((theta_results[i] - vG_WCR[i])/ (vG_WCS[i] - vG_WCR[i]))

#--------------------------------------------------
#dat file 
#--------------------------------------------------
#Stress - strain------
np.savetxt(f"../results/stress_strain.dat",
    np.column_stack((sigma_eff, -strain)),
    header="sigma_eff_kPa strain",fmt="%.8e")

#Pressure - void ratio------
np.savetxt(f"../results/pressure_void_simulation.dat",
    np.column_stack((Pressure_results[i],e_results[i])),
    header="pore_pressure_kPa void_ratio",fmt="%.8e")

#Pressure - strain------------
np.savetxt(f"../results/pressure_strain_simulation.dat",
    np.column_stack((Pressure_results[i],-strain)),
    header="pore_pressure_kPa strain",fmt="%.8e")

#Void ratio - time----------
np.savetxt(f"../results/void_time_simulation.dat",
    np.column_stack((simulation_time,e_results[i])),
    header="time_day void_ratio",fmt="%.8e")

# Pressure - time----------
np.savetxt(f"../results/pressure_time_simulation.dat",
    np.column_stack((simulation_time,Pressure_results[i])),
    header="time_day pore_pressure_kPa",fmt="%.8e")

# Soil shrinkage characteristic curve - Se from theta-------
np.savetxt(f"../results/sscc_vG.dat",
    np.column_stack((Se_theta,e_results[i])),
    header="Se void_ratio",fmt="%.8e")
# Moisture ratio - void ratio -------
np.savetxt(f"../results/sscc_MR_simulation.dat",
    np.column_stack((MR[i],e_results[i])),
    header="moisture_ratio void_ratio",fmt="%.8e")


#--------------------------------------------------
#Plotting 
#--------------------------------------------------
plt.figure()
plt.plot(sigma_eff,-strain,".",color="black")
plt.title(f"Stress strain")
plt.xlabel(rf"Effective stress $\sigma'$[kPa] ")
plt.ylabel(rf"Strain $\varepsilon$[-]")
plt.xscale("log")
textstr = (
    rf"$\kappa = {Compression_vector.Kappa[i]:.3f}$" "\n"
    rf"$\lambda = {Compression_vector.Lambda[i]:.3f}$" "\n"
    rf"$\mu = {Compression_vector.Mu[i]:.3f}$" "\n"
    rf"$\mathrm{{OCR}} = {Compression_vector.OCR[i]:.3f}$" "\n"
    rf"$C = {Compression_vector.c[i]} $" "\n"
    rf"$\phi = {Compression_vector.phi[i]}^\circ$"
)

plt.text(0.05,0.05,textstr,transform=plt.gca().transAxes,fontsize=11,
    verticalalignment="bottom",bbox=dict(boxstyle="round",facecolor="white",alpha=0.8))
plt.grid(True, which="both", linestyle="--", alpha=0.5)
plt.savefig('../results/figure/stress_strain.png')

plt.figure()
plt.plot(Pressure_results[i], e_results[i],color="red")
plt.title(f"Stress strain at i ={i}")
plt.xlabel(rf"Pore pressure $\psi$ [kPa]")
plt.ylabel(rf"Void ratio $e$ [-]")
#plt.xscale("log")
plt.grid(True, which="both", linestyle="--", alpha=0.5)
plt.savefig('../results/figure/pressure_void.png')

plt.figure()
plt.plot(Pressure_results[i], -strain,color="red")
#plt.title(f"Stress strain at i ={i}")
plt.xlabel(rf"Pore pressure $\psi$ [kPa]")
plt.ylabel(rf"Strain $\varepsilon$ [-]")
#plt.xscale("log")
plt.grid(True, which="both", linestyle="--", alpha=0.5)
plt.savefig('../results/figure/pressure_strain.png')

#Void_time------
plt.figure()
plt.plot(simulation_time,e_results[i])
plt.title(f"Pressure at i ={i}")
plt.xlabel("Time [day]")
plt.ylabel(rf"Void ratio $e$")
plt.savefig('../results/figure/void_time.png')

#Pressure_time------
plt.figure()
plt.plot(simulation_time,Pressure_results[i])
plt.title(f"Pressure at i ={i}")
plt.xlabel("Time [day]")
plt.ylabel("Pore pressure [kPa]")
plt.savefig('../results/figure/pressure_time.png')

#Soil shrinkage characteristic curve------
plt.figure()
plt.plot(MR[i],e_results[i],label="Simulation")
#plt.title(f"Soil shrinkage characteristic at i ={i}")
plt.xlabel(rf"Moisture ratio $\vartheta$")
plt.ylabel(rf"Void ratio $e$")
plt.grid(True, which="both", linestyle="--", alpha=0.5)
plt.savefig('../results/figure/sscc_MR.png')

#Groundwater level------
#gwl = []
#gwl_check = 0
#for i in range(Pressure_results.shape[1]):  # Loop over time
#    for j in range(Pressure_results.shape[0] - 1):  # Loop over depth
#        if Pressure_results[j, i] > 0:
#            gwl.append(colum_node_height[j, i] + Soil_dimension.Depth_bottom)
#            gwl_check = 1
#            break 
#
#if gwl_check == 1:
#    np.savetxt(f"{dat_dir}/gwl.dat",
#        np.column_stack((simulation_time,gwl)),
#        header="Time[day] Groundwater_level[m]",fmt="%.8e")
