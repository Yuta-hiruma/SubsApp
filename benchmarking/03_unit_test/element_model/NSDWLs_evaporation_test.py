# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 11:25:34 2024

@author: tomdegast
"""

import matplotlib.pyplot as plt
import timeit as timeit
import pandas as pd
import math as math
import scipy as sci
import numpy as np
from scipy import special
from NSDWLs_flux import *
from NSDWLs_calc import *
from NSDWLs_potential import *
from NSDWLs_equations import *
from NSDWLs_states import *
from NSDWLs_errors import *
from build_domain import *
from Input_parameters_definitiions import *
" ## local files##"
# from main_program import *


" ## global files ##"

" #### input parameters "

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
Compression_vector = Compression_vector()
Compression_vector.ccott = 30.
Compression_vector.OCR = 1.2
Compression_vector.Lambda = 0.114
Compression_vector.Kappa = 0.078
Compression_vector.Mu = 0.045
Compression_vector.Tau_ref = 1
Compression_vector.unloadreload = 0.9

Compression_vector.K0 = 0.1
Compression_vector.poisson = 0.1

Water_retention_curve_node = Water_retention_curve_vectors()
Water_retention_curve_node.vG_alpha = 0.11
Water_retention_curve_node.vG_n = 1.23
Water_retention_curve_node.vG_m = (1-1/1.23)
Water_retention_curve_node.vG_a = 50000
Water_retention_curve_node.vG_mod = 1

PP = 400
e0 = 2.87634
h0 = 0.02
S_e = 1

j = 0
for i in range(len(simulation_time)):
    if time_incr[j] < simulation_time[i]:
        time_incr_step[j] = simulation_time[i]
        j = j+1
        if j == len(time_incr):
            break


dsigmaeff_dt = 0
sigma_tot = 12*(h0/2)
u = 10*(h0/2)
sigma_eff = sigma_tot-u


def equations(simulation_time, in_in):
    global load_inc
    global load
    global time_pressure
    global pressure
    global time_incr
    global Compression_vector
    global dload_increment
    global time_incr_step
    global time_step
    global h0
    global e0

    " get "
    PP_temp = in_in.copy
    eps = in_in[0]
    PP = in_in[1]
    sigma_eff = in_in[2]
    u = in_in[3]
    simga_tot = in_in[4]
    S_e = in_in[5]
    e = in_in[6]


    print("simulation_time", round(simulation_time, 2))


    " calculate void ratio"
    #e = (1+e0)*(eps+1)-1
    # e = (1+e0)*(1+eps)-1
    if S_e > 1:
        S_e = 1

    theta = S_e*e


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

    # if simulation_time >= time_incr[t_interval] + time_step:
        # sigma_eff = load[t_interval]

    " get pressure boundary conditions"
    t_interval = len(time_pressure)-1

    for i in range(len(time_pressure)-1):
        # print(i)
        if simulation_time < time_pressure[i+1]:
            t_interval = t_interval-1

    u_t = -pressure[t_interval]

    if simulation_time <= time_pressure[t_interval]:
        u_t = -pressure[t_interval]

    dsimga_tot_dt = -dsigmaeff_dt
    #dsimga_tot_dt = dsigmaeff_dt

    "Calculate permeability based on void ratio"
    # S_e = 1

    K = 10**(0.4*(S_e*e)-4)
    K = 1

    Water_retention_curve_vector = 1

    " set change in pore pressure"
    du_dt = K*d2f_dz2_non_uniform_loc(u, u_t, u_t, h0)
    "check compression or extension"
    compression = 1

    if dsigmaeff_dt < 0:
        compression = 0

    "check if saturated"

    saturated = 1

    # if S_e <1 : saturated =0

    " add total boundary load to water pressure"
    du_dt = du_dt - dsimga_tot_dt  # *S_e #(-1/K)*compression*saturated*deps_dt

    "calculate saturated effective stress"

    #dsigmaeff_dt = (dsimga_tot_dt + du_dt)*S_e
    dsigmaeff_dt = (dsimga_tot_dt - du_dt)*S_e

    #print("dsigmaeff_dt", round(dsigmaeff_dt, 2))

    [deps_dt, dPP_dt] = deSSc_dt_loc_h(Compression_vector, -sigma_eff, dsigmaeff_dt,e, Water_retention_curve_vector, PP)

    de_dt = vertical_strain_to_void_ratio(Compression_vector,e,deps_dt)


    dS_e_dpsi = derrivative_van_genuchten_S_e_loc(Water_retention_curve_node, theta, e, psi_m)

    # dSat_dt = dtheta_dt/e - unsaturated*de_dt/e
    dSat_dt = -dS_e_dpsi*du_dt

    return [deps_dt, dPP_dt, dsigmaeff_dt, du_dt, dsimga_tot_dt, dSat_dt, de_dt]


input_in = [0, PP, sigma_eff, u, sigma_tot, S_e, e0]

sol = solve_ivp(equations, [simulation_time[0], simulation_time[-1]], input_in, method='RK45',
                t_eval=simulation_time, max_step=max_time_step, rtol=rtol)  # , atol = atol)

with open('output_evaporation/void_time.dat', 'w') as dis:
    dis.write(f"{'':<2}{'time[d]':<15}  {'void ratio':<15}\n")
    for i in range(len(simulation_time)):
        dis.write(f"{'':<2}{simulation_time[i]:.4f}  {sol.y[6,i]:.4f}\n")


#plt.figure()
#plt.scatter(suction_measured, e_measured, marker='.')
#plt.ylabel("void ratio")
#plt.xlabel("pressure")
#plt.scatter(-sol.y[3, :], (1+e0)*(sol.y[0, :]+1)-1)
#plt.xlim(0.1, 10000)
#plt.ylim(0, 4)
#plt.text(1000, 3.75, f"Lambda: {Compression_vector.Lambda:.4f}")
#plt.text(1000, 3.65, f"Kappa: {Compression_vector.Kappa:.4f}")
#plt.text(1000, 3.55, f"Mu: {Compression_vector.Mu:.4f}")
#plt.text(1000, 3.45, f"OCR: {PP:.2f}")
#plt.text(1000, 3.35, f"ccott: {Compression_vector.ccott:.2f}")
#plt.xscale('log')
#plt.legend()
#plt.savefig('output_evaporation/linear_strain.png')



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


