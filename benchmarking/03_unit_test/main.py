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
import numpy as np
import scipy as sci
from scipy.integrate import solve_ivp
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
from Input_parameters_definitiions import *
#Reading input data ----------
import sys
import os


#Import laboratory test data -----------------------

file_path = 'test/hysterese_OV_plot.csv'
try:
    df = pd.read_csv(file_path, sep=';')
except FileNotFoundError:
    print("file doesn't exist")

time_values = df['time'].to_numpy()
suction_measured = df['suction'].to_numpy()
e_measured = df['void ratio'].to_numpy()

dsuction = np.zeros(len(e_measured)-1)


#Simulaiton time ------------------------------------------------
time_pressure = time_values
pressure = suction_measured

time_incr = [0]
time_incr_step = np.zeros(len(time_incr))
time_max = max(time_values)
#time_max = 25

simulation_time = np.linspace(0, time_max, 2500)
time_step = simulation_time[1]
max_time_step = time_step
rtol = 1e-4  # relative tolerance

#No boundary condition in this element test-----


#Initial condition --------------------------------------------
e0 = 2.87634
h0 = 0.02
S_e = 1
material = "Clay"

dsigmaeff_dt = 0
sigma_tot = 12*(h0/2)
u = 10*(h0/2)
sigma_eff = sigma_tot-u

Soil_stresses.psi_ini = u
Soil_states.theta = e0
Soil_states.eps_ini = 0.0

Soil_dimension = Soil_dimension()
Soil_dimension.node_distance = h0/2
Bctop_in = []
Bctop_in.pressure = []
#Material parameter --------------------------------------------
gamma_vector = gamma_vectors()
gamma_vector.s = Soil_parameters[material]["Material_parameters"]["gamma_s"]
gamma_vector.w = Soil_parameters[material]["Material_parameters"]["gamma_w"]

Permeability_curve_Perm_vector = Permeability_curve_Perm_vector()
Permeability_curve_Perm_vector.K_tot = Soil_parameters[material]["Permeability_curve_Perm"]["K_tot"]
Permeability_curve_Perm_vector.K_min = Soil_parameters[material]["Permeability_curve_Perm"]["K_min"]
Permeability_curve_Perm_vector.tau = Soil_parameters[material]["Permeability_curve_Perm"]["tau"]

Compression_vector = Compression_vector()
Compression_vector.Kappa =Soil_parameters[material]["Compression"]["kappa"]
Compression_vector.Lambda = Soil_parameters[material]["Compression"]["lambda"]
Compression_vector.Mu = Soil_parameters[material]["Compression"]["mu"]
Compression_vector.Tau_ref = Soil_parameters[material]["Compression"]["Tau_ref"]
Compression_vector.Pg = Soil_parameters[material]["Compression"]["Pg"]
Compression_vector.OCR = Soil_parameters[material]["Compression"]["OCR"]
Compression_vector.c = Soil_parameters[material]["Compression"]["c"]
Compression_vector.phi = Soil_parameters[material]["Compression"]["phi"]
Compression_vector.ccott = Compression_vector.c/math.tan(math.radians(Compression_vector.phi))
#Horizontal volume change ------------
Compression_vector.K0 = 0.1
Compression_vector.poisson = 0.1

Water_retention_curve_vector = Water_retention_curve_vectors()
Water_retention_curve_vector.vG_WCR = Soil_parameters[material]["Water_retention_curve"]["vG_WCR"]
Water_retention_curve_vector.vG_WCS = Soil_parameters[material]["Water_retention_curve"]["vG_WCS"]
Water_retention_curve_vector.vG_a = Soil_parameters[material]["Water_retention_curve"]["vG_a"]
Water_retention_curve_vector.vG_alpha = Soil_parameters[material]["Water_retention_curve"]["vG_alpha"]
Water_retention_curve_vector.vG_m = Soil_parameters[material]["Water_retention_curve"]["vG_m"]
Water_retention_curve_vector.vG_mod = Soil_parameters[material]["Water_retention_curve"]["vG_mod"]
Water_retention_curve_vector.vG_n = Soil_parameters[material]["Water_retention_curve"]["vG_n"]

#loading ------------------------------------------------------
load_increments = [0]
load = np.cumsum(load_increments)
dload_increment = load_increments/simulation_time[1]


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
#Initialize tracker -----------------------------------
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


    #No boundary condition in this element test-----
                                                       
 
    #Update soil_dimension------------
    #Soil_dimension.node_vector = get_node_coordinates(Soil_dimension)
    #Soil_dimension.node_distance = get_node_distance(Soil_dimension)  
    #Soil_dimension.node_distance = (1-((Soil_states.e_ini-e)/(1+Soil_states.e_ini)))*Soil_dimension.node_initial
 
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
    K= Soil_calc.K/Soil_dimension.node_distance

    #set differentials ------------------------
    #First set: pressure ----------
    d2psi_dz2 = d2f_dz2_non_uniform(Pressure,Bctop_in.pressure,Bcbottom_in.pressure,Soil_dimension) #Pressure*0+np.sin(simulation_time)
    dpsi_dt = K*d2psi_dz2
    dS_e_dpsi,dtheta_dpsi  = derrivative_van_genuchten_S_e(Water_retention_curve_node, theta,e, psi_m)
    # second set: theta ----------
    dtheta_dt = dpsi_dt*dtheta_dpsi  
    # third set: void ratio -------
    de_dt,dn_dt = deSSc_dt(Compression_vector, Soil_states, S_e, sigma_eff, dsigma_eff_dt,dpsi_dt,Water_retention_curve_vector, Soil_states.Preconsolidation_pressure)
    de_dt = vertical_strain_to_void_ratio(Compression_vector,e,deps_dt)


    compression = (de_dt < 0).astype(int)    
    extension = (de_dt > 0).astype(int)  
    
    positive_pressure = (dpsi_dt > 0).astype(int)
    negative_pressure = (dpsi_dt < 0).astype(int)

    dpsi_dt = update_dpsi_dt(dpsi_dt, K, compression, saturated, dn_dt)
    dtheta_dt, dSat_dt = update_dtheta_and_dSat(dtheta_dt,positive_pressure, extension, saturated, dn_dt,n, theta, unsaturated)
    dtheta_dt =  gather(*(saturated_theta(S_e,dtheta_dt,dn_dt,i) for i in range(len(S_e))))#.compute()

    out_out = np.concatenate([dpsi_dt,dtheta_dt,de_dt],axis=0)
    return [out_out]



#Set input --------- # Calculate strain instead of void ratio#
# First set is pressure #
input_in = Soil_stresses.psi_ini
# Second set is water content #
input_in = np.append(input_in,Soil_states.theta)
# Third set is strain #
input_in = np.append(input_in,Soil_states.eps_ini) 


# Solveing differential equation ------
start = time.perf_counter() # measure the performance
sol = solve_ivp(equations, [simulation_time[0], simulation_time[-1]],input_in,method='RK23',t_eval=simulation_time,max_step = max_time_step, rtol = rtol ,args=(sigma_tracker,))#, atol = atol)
end = time.perf_counter()# measure the performance

print("Elapsed  = {}s".format((end - start)))
print("Success:", sol.success)
print("Message:", sol.message)
print("Final time reached:", sol.t[-1])

#Post process---------
nr_elements = int(np.size(input_in)/3)
# First set is pressure -----
Pressure_results = sol.y[0:nr_elements]
# Second set is water content-----
theta_results = sol.y[nr_elements:2*nr_elements]
# Third set is void ratio-----
e_results = sol.y[2*nr_elements:3*nr_elements]
#-------------------------------------------------------------
n_results = e_results/(1+e_results)
#Se_results = theta_results/n_results
Se_results = theta_results/e_results

#Convert void ratio to surface displacement -----
colum_heigth=np.zeros(len(e_results[0,:]))
node_heigth=0*(e_results)
for i in range (len(e_results[0,:])):
    colum_heigth[i] = np.sum((1-((Soil_states.e_ini-e_results[:,i])/(1+Soil_states.e_ini)))*Soil_dimension.node_initial)
    node_heigth[:,i] = ((1-((Soil_states.e_ini-e_results[:,i])/(1+Soil_states.e_ini)))*Soil_dimension.node_initial)
colum_node_heigth = node_heigth*0
colum_node_heigth[-1,:] = node_heigth[-1,:]

for i in range(len(e_results[:,0])-2,-1,-1):
    colum_node_heigth[i,:] = node_heigth[i,:]+colum_node_heigth[i+1,:]

#ground water level-------------------------------------------
gwl = []

for i in range(Pressure_results.shape[1]):  # Loop over columns (time)
    for j in range(Pressure_results.shape[0] - 1):  # Loop over rows (depth)
        if Pressure_results[j, i] > 0:
            gwl.append(colum_node_heigth[j, i] + Soil_dimension.Depth_bottom)
            break 


with open('output_evaporation/void_time.dat', 'w') as dis:
    dis.write(f"{'':<2}{'time[d]':<15}  {'void ratio':<15}\n")
    for i in range(len(simulation_time)):
        dis.write(f"{'':<2}{simulation_time[i]:.4f}  {sol.y[6,i]:.4f}\n")

plt.figure()
plt.scatter(suction_measured, e_measured, marker='.', label='Laboratory test')
plt.ylabel("void ratio [-]")
plt.xlabel("Suction [kPa]")
plt.scatter(-sol.y[3, :], sol.y[6,:],label='Simulation')
plt.xlim(0.1, 10000)
plt.ylim(0, 4)
#plt.text(1000, 3.75, f"Lambda: {Compression_vector.Lambda:.3f}")
#plt.text(1000, 3.63, f"Kappa: {Compression_vector.Kappa:.3f}")
#plt.text(1000, 3.51, f"Mu: {Compression_vector.Mu:.3f}")
#plt.text(1000, 3.39, f"OCR: {Compression_vector.OCR:.1f}")
#plt.text(1000, 3.27, f"ccott: {Compression_vector.ccott:.0f}")
#plt.text(1000, 3.15, f"poisson ratio: {Compression_vector.poisson:.1f}")
#plt.text(1000, 3.03, f"K0: {Compression_vector.K0:.1f}")
plt.xscale('log')
plt.legend()
plt.savefig('output_evaporation/linear_strain_3D.png')
plt.figure()


#Output ------------------------------------------
with open('../results/displacement.dat', 'w') as dis:
    dis.write(f"{'':<2}{'time[days]':<15}  {'displacement [mm]':<15}\n")
    for i in range(len(date_list)):
        dis.write(f"{'':<2}{date_list[i].strftime('%Y-%m-%d'):<15}  {(colum_heigth[i]+Soil_dimension.Depth_bottom)*1000:<15}\n")

with open('../results/Pressure_top.dat', 'w') as dis:
    dis.write(f"{'':<2}{'time[days]':<15} {'Pressure top':<15} \n")
    for i in range(len(date_list)):
        dis.write(f"{'':<2}{date_list[i].strftime('%Y-%m-%d'):<15} {Pressure_results[0,i]:<15} \n")

with open('../results/Pressure_bottom.dat', 'w') as dis:
    dis.write(f"{'':<2}{'time[days]':<15} {'Pressure top':<15} \n")
    for i in range(len(date_list)):
        dis.write(f"{'':<2}{date_list[i].strftime('%Y-%m-%d'):<15} {Pressure_results[-1,i]:<15} \n")

with open('../results/Se.dat', 'w') as dis:
    dis.write(f"{'':<2}{'time[days]':<15} {'depth':<15} {'Se':<15}\n")
    for i in range(len(date_list)):
        dis.write(f"{'':<2}{date_list[i].strftime('%Y-%m-%d'):<15} {colum_node_heigth[1,i]+Soil_dimension.Depth_bottom:<15} {theta_results[1,i]*(1+e_results[1,i])/e_results[1,i]:<15}\n")

with open('../results/theta.dat', 'w') as dis:
    dis.write(f"{'':<2}{'time[days]':<15} {'depth':<15}  {'theta':<15}\n")
    for i in range(len(date_list)):
        dis.write(f"{'':<2}{date_list[i].strftime('%Y-%m-%d'):<15} {colum_node_heigth[1,i]+Soil_dimension.Depth_bottom:<15} {theta_results[1,i]:<15}\n")

with open('../results/n.dat', 'w') as dis:
    dis.write(f"{'':<2}{'time[days]':<15}  {'depth':<15} {'n':<15}\n")
    for i in range(len(date_list)):
        dis.write(f"{'':<2}{date_list[i].strftime('%Y-%m-%d'):<15} {colum_node_heigth[1,i]+Soil_dimension.Depth_bottom:<15} {e_results[1,i]/(1+e_results[1,i]):<15}\n")

#Plot results----------------------------------
plt.figure()
plt.plot(date_list,(colum_heigth+Soil_dimension.Depth_bottom)*1000)
plt.xlabel('time [days]')
plt.ylabel('displacement [mm]')
plt.title('displacement')
plt.xticks(rotation=45)
plt.subplots_adjust(left=0.2, bottom=0.2)
plt.savefig('../results/figure/displacement.png')

#--------------
abs_max = np.max(np.abs(Pressure_results))
#abs_min = np.min(np.abs(Pressure_results))
plt.figure()
plt.pcolormesh(date_list, Soil_dimension.node_vector, Pressure_results,
               cmap='bwr', vmin=-abs_max, vmax=abs_max)
for i in range(len(Soil_dimension.Depth_o_soils) - 1):
    plt.axhline(y=Soil_dimension.Depth_o_soils[i], color='k', linestyle='--')
plt.title("Pressure")
plt.xlabel("Time [days]")
plt.ylabel("Depth [m]")
plt.xticks(rotation=45)
plt.subplots_adjust(left=0.2, right=0.8, bottom=0.2)
cax = plt.axes((0.85, 0.1, 0.075, 0.8))
plt.colorbar(cax=cax)
plt.title("kPa")
plt.savefig('../results/figure/pressure.png')
#--------------
plt.figure()
#plt.pcolormesh(sol.t,Soil_dimension.node_vector,theta_results, cmap='RdBu', vmin=np.min(theta_results), vmax=np.max(theta_results))
plt.pcolormesh(date_list,Soil_dimension.node_vector,theta_results, cmap='Blues', vmin=np.min(theta_results), vmax=np.max(theta_results))
for i in range(len(Soil_dimension.Depth_o_soils)-1):
    plt.axhline(y=Soil_dimension.Depth_o_soils[i], color='k', linestyle='--')
plt.title("theta")
plt.xticks(rotation=45)
plt.subplots_adjust(left=0.2, bottom=0.2) 
plt.xlabel("time [days]")
plt.ylabel("depth [m]")
plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
cax = plt.axes((0.85,0.1,0.075,0.8))
plt.colorbar(cax=cax)
plt.title("[-]")
plt.savefig('../results/figure/theta.png')

#--------------
plt.figure()
vG_WCR = Water_retention_curve_vector.vG_WCR[:, np.newaxis]  
vG_WCS = Water_retention_curve_vector.vG_WCS[:, np.newaxis]
norm = mcolors.Normalize(vmin=np.min((theta_results-vG_WCR)/(vG_WCS-vG_WCR)), vmax=np.max((theta_results-vG_WCR)/(vG_WCS-vG_WCR)))
plt.pcolormesh(date_list,Soil_dimension.node_vector,(theta_results-vG_WCR)/(vG_WCS-vG_WCR), cmap='Blues', norm=norm)
for i in range(len(Soil_dimension.Depth_o_soils)-1):
    plt.axhline(y=Soil_dimension.Depth_o_soils[i], color='k', linestyle='--')
plt.title("Se")
plt.xlabel("time [days]")
plt.ylabel("depth [m]")
plt.xticks(rotation=45)
plt.subplots_adjust(left=0.2, bottom=0.5) 
plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
cax = plt.axes((0.85,0.1,0.075,0.8))
plt.colorbar(cax=cax)
plt.title("[-]")
plt.savefig('../results/figure/Se_vG.png')

#--------------
plt.figure()
norm = mcolors.Normalize(vmin=np.min(e_results/(1+e_results)), vmax=np.max(e_results/(1+e_results)))
plt.pcolormesh(date_list, Soil_dimension.node_vector, e_results/(1+e_results), cmap='Blues', norm=norm)
for i in range(len(Soil_dimension.Depth_o_soils)-1):
    plt.axhline(y=Soil_dimension.Depth_o_soils[i], color='k', linestyle='--')
plt.title("Porosity")
plt.xlabel("time [days]")
plt.ylabel("depth [m]")
plt.xticks(rotation=45)
plt.subplots_adjust(left=0.2, bottom=0.2) 
plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
cax = plt.axes((0.85,0.1,0.075,0.8))
plt.colorbar(cax=cax)
plt.title("[-]")
plt.savefig('../results/figure/porosity.png')
