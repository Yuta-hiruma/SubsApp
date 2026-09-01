# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 11:25:34 2024

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

for i in range(len(suction_measured)-1):
    dsuction[i] = (suction_measured[i]-suction_measured[i+1]) / \
        (time_values[i]-time_values[i+1])

load_increments = [0]
load = np.cumsum(load_increments)

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


dload_increment = load_increments/simulation_time[1]

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
Soil_states.theta = e0/(1+e0) 
Soil_states.eps_ini = 0.0

Soil_dimension = Soil_dimension()
Soil_dimension.node_distance = float(h0/2)

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
Compression_vector.OCR = Soil_parameters[material]["Compression"]["OCR"]
#Compression_vector.c = Soil_parameters[material]["Compression"]["c"]
#Compression_vector.phi = Soil_parameters[material]["Compression"]["phi"]
#Compression_vector.ccott = Compression_vector.c/math.tan(math.radians(Compression_vector.phi))
Compression_vector.ccott = Soil_parameters[material]["Compression"]["ccott"]
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


j = 0
for i in range(len(simulation_time)):
    if time_incr[j] < simulation_time[i]:
        time_incr_step[j] = simulation_time[i]
        j = j+1
        if j == len(time_incr):
            break


def equations(simulation_time, in_in):
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
        
        Returns
        -------
        ndarray
            Time derivatives of the state variables (psi, theta, and e) for all
            soil elements.
    """    

    global load_inc
    global load
    global time_pressure
    global pressure
    global time_incr
    global gamma_vector
    global Compression_vector
    global Water_retention_curve_vector
    global Permeability_cruve_Perm_vector
    global Soil_dimension
    global dload_increment
    global time_incr_step
    global time_step
    global h0
    global e0

    u = in_in[0]
    Soil_states.theta = in_in[1]
    Soil_states.e = in_in[2]

    print("simulation_time", round(simulation_time, 2))


    S_e = np.clip(update_S_e(Water_retention_curve_vector, Soil_states.theta ), 0.0, 1.0)
    saturated = (S_e >= 1).astype(int)
    unsaturated = (S_e < 1).astype(int)  

    " get suction forom u"
    psi_m = np.where(u <= 0, u, 0)

    " get load boundary conditions"
    t_interval = len(time_incr)-1

    for i in range(len(time_incr)):
        # print(i)
        if simulation_time < time_incr[i]:
            t_interval = t_interval-1

    dsigmaeff_dt = 0

    if simulation_time <= time_incr[t_interval] + time_step:
        dsigmaeff_dt = -dload_increment[t_interval]
    #get pressure boundary conditions--
    t_interval = len(time_pressure)-1

    for i in range(len(time_pressure)-1):
        if simulation_time < time_pressure[i+1]:
            t_interval = t_interval-1

    u_t = -pressure[t_interval]

    if simulation_time <= time_pressure[t_interval]:
        u_t = -pressure[t_interval]

    #Update_soil_stresses------------
    #sigma_total = update_sigma_total (Soil_dimension.node_distance,gamma_vector,Soil_states,Water_retention_curve_vector)
    #sigma_eff = np.abs(sigma_total - Pressure*S_e)

    dsimga_tot_dt = -dsigmaeff_dt
    #dsimga_tot_dt = dsigmaeff_dt

    #Calculate permeability----------
    #K = 10**(0.4*(S_e*e)-4)
    #K = 1
    K = update_K_Kmin (Soil_states.theta,Soil_states.e,S_e,Permeability_curve_Perm_vector,Water_retention_curve_vector)


    " set change in pore pressure"
    du_dt = K*d2f_dz2_non_uniform_loc(u, u_t, u_t, h0)

    compression = 1
    if dsigmaeff_dt < 0:
        compression = 0
    du_dt = du_dt - dsimga_tot_dt  # *S_e #(-1/K)*compression*saturated*deps_dt

    #dsigmaeff_dt = (dsimga_tot_dt + du_dt)*S_e
    dsigmaeff_dt = (dsimga_tot_dt - du_dt)*S_e

    deps_dt = deSSc_dt_loc(Compression_vector, -sigma_eff, dsigmaeff_dt, Soil_states.e, Water_retention_curve_vector)
    de_dt = vertical_strain_to_void_ratio(Compression_vector, Soil_states.e, deps_dt)
    #de_dt = deps_dt

    dS_e_du, dtheta_du = derrivative_van_genuchten_S_e(Water_retention_curve_vector, Soil_states.theta, Soil_states.e, psi_m)
    dtheta_dt = -du_dt*dtheta_du 


    return [du_dt, dtheta_dt, de_dt]

#input_in = [0, PP, sigma_eff, u, sigma_tot, S_e, e0]
input_in = [u, Soil_states.theta, e0]

# Solveing differential equation ------
start = time.perf_counter() # measure the performance
sol = solve_ivp(equations, [simulation_time[0], simulation_time[-1]],input_in,method='RK23',t_eval=simulation_time,max_step = max_time_step, rtol = rtol)#, atol = atol)
end = time.perf_counter()# measure the performance

with open('output_evaporation/void_time.dat', 'w') as dis:
    dis.write(f"{'':<2}{'time[d]':<15}  {'void ratio':<15}\n")
    for i in range(len(simulation_time)):
        dis.write(f"{'':<2}{simulation_time[i]:.4f}  {sol.y[2,i]:.4f}\n")

#Plot void results----------------------------------
#plt.figure()
#fig, ax = plt.subplots(figsize=(8, 5))
#ax.plot(simulation_time,sol.y[2],label="Simulated",linewidth=1,color="blue" )
#ax.scatter(time_values,e_measured,label="Measured",color="orange",linewidth=1,zorder=3,)
#ax.set(xlabel="Time [days]",ylabel="Void ratio [-]",title="Void ratio over time",)
#ax.grid(True, alpha=0.3)
#ax.legend()
#fig.autofmt_xdate(rotation=45)
#fig.tight_layout()
#fig.savefig('output_evaporation/displacement.png',dpi=300, bbox_inches="tight",)



plt.figure()
plt.scatter(suction_measured, e_measured, marker='.', label='Laboratory test')
plt.ylabel("void ratio [-]")
plt.xlabel("Suction [kPa]")
plt.scatter(-sol.y[0, :], sol.y[2,:],label='Simulation')
plt.xlim(0.1, 10000)
plt.ylim(0, 4)
plt.xscale('log')
plt.legend()
plt.savefig('output_evaporation/linear_strain_3D.png')
plt.figure()


plt.figure()
#plt.scatter(suction_measured, e_measured, marker='.', label='Laboratory test')
plt.xlabel("Se [-]")
plt.ylabel("Void ratio [-]")
plt.scatter(sol.y[1, :], sol.y[2,:],label='Simulation')
#plt.xlim(0.1, 10000)
#plt.ylim(0, 4)
#plt.xscale('log')
plt.legend()
plt.savefig('output_evaporation/sscc.png')
plt.figure()

