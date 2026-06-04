# -*- coding: utf-8 -*-
"""
Created on Fri Jul 23 14:27:28 2021

@author: tomdegast

This model links the change in effective stress due to transient flow to a consitutitve modelincluding creep (Modelled after Soft Soil creep). 
The permeability is changed based on watercontent and void ratio. 


"""
"Numba"
from numba import int32, float32
from numba import jit
import time

"The following function libraries are used in the programme and are local"
" ## local files##" 

#from main_program import *
from build_domain import *
from NSDWLs_errors import *
from NSDWLs_states import *
from NSDWLs_equations import *
from NSDWLs_potential import *


"The following function libraries are used in the programme and are global"

" ## global files ##"
from scipy import special
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as ticker 
import matplotlib.colors as mcolors
import numpy as np
import scipy as sci
import math as math
import pandas as pd
import timeit as timeit
from datetime import datetime, timedelta

#------
import sys
import os
current_dir = os.getcwd()
parent_dir = os.path.abspath("../data")
sys.path.insert(0, parent_dir)
from Input_soil_material_parameters import *
from input import *
sys.path.pop(0)
os.chdir(current_dir)
#------

#input file-----------------------------
infile_tp = "input/top_pressure/top_pressure_swp.dat"
infile_bp = "input/bottom_pressure/bottom_pressure_swp.dat"
infile_wc_ini = "input/wc_ini/wc_ini_swp.dat"
infile_p_ini = "input/pressure_ini/pressure_ini_swp.dat"
df_tp = pd.read_csv(infile_tp, sep=r"\s+", skiprows=1, names=["date", "pressure"])
df_bp = pd.read_csv(infile_bp, sep=r"\s+", skiprows=1, names=["date", "pressure"])
df_wc_ini = pd.read_csv(infile_wc_ini, sep=r"\s+", skiprows=1, names=["depth", "wcontent"])
df_pressure_ini = pd.read_csv(infile_p_ini, sep=r"\s+", skiprows=1, names=["depth", "pressure"])
#---------------------------------------
#simulaiton period----------------------
TSTART = tstart
TEND = tend

" #### input parameters "
tic = timeit.default_timer()

" Initiation of he soil domain, give top, bottom and precision"
" precision directly influences calculation time. If the precision is small, the relative error required for accuracy will increase calulation time"
Soil_dimension = Soil_domain()
Soil_dimension.Depth_top = 0.00                     # m size
Soil_dimension.Depth_bottom = depths[-1]            # m size
Soil_dimension.precision = 0.1                      # m node size

" Layers are added to the domain."
" Soil properties are collected in the library Input_parameters_definitiions "

Soil_dimension.Depth_o_soils = depths                # bottom of soil layer      
Soil_dimension.soils = types                                      # soil types used

""

Soil_dimension.Number_o_soils = len(Soil_dimension.Depth_o_soils)                                                 # number of soils 
Soil_dimension.domainsize = Soil_dimension.Depth_top-Soil_dimension.Depth_bottom

Soil_states = Soil_states()
Soil_stresses = Soil_stresses()

Soil_calc = Soil_calc()
swp = swp()

" ## Set calculation time, simulation steps to be reported, maximum time step and relative tolerance ##"
start_date = datetime.strptime(TSTART, "%Y-%m-%d")
end_date = datetime.strptime(TEND, "%Y-%m-%d")
minimum_accuracy_time_step = 1               #days  
max_time_step = 1.0  
rtol =1e-7 #relative tolerance
#atol =1e-4 #absolute tolerance

end_time = (end_date - start_date).days
simssteps = (end_date - start_date).days
date_list = [start_date + timedelta(days=i) for i in range((end_date - start_date).days + 1)]
simulation_time =  np.linspace (0,end_time,simssteps+1)



"## Set boundary conditions"

"Top boundary conditions"
Bctop = Bc_top()
#Initial conditoin from SWAP----------
Bctop.time = df_tp["date"].tolist()
Bctop.pressure = df_tp["pressure"].tolist()   
#-------------------------------------
Bctop.trans = np.ones(len(Bctop.pressure))
Bctop.load_increments = [0]
Bctop.load = np.cumsum(Bctop.load_increments)
Bctop.time_incr = [0]


"Bottom boundary conditions"        
Bctop.devaporation = np.zeros(len(Bctop.pressure))
Bcbottom = Bc_bottom()
#Bcbottom.time = [0]
#Bcbottom.pressure = [49]
Bcbottom.time = df_bp["date"].tolist()
Bcbottom.pressure = df_bp["pressure"].tolist()     # kPa

"## check for errors ##"

checkerror(Soil_dimension,Soil_parameters)

" Build the model domain and fill vectors with properties set initial conditions" 
#Initial conditoin from SWAP----------
swp.depth = df_wc_ini ["depth"].tolist()
swp.wcontent = df_wc_ini ["wcontent"].tolist()
swp.pressure = df_pressure_ini ["pressure"].tolist() #kPa
#-------------------------------------
(Soil_states,Soil_states.theta_ini,Soil_stresses.psi_ini) = set_ini(Soil_dimension,Soil_parameters,Soil_states,Bctop,swp)
Soil_states.e_ini = Soil_states.e
Soil_states.n_ini = Soil_states.n
Soil_states.theta = Soil_states.theta_ini
Soil_dimension.internode_vector = modeldomain (Soil_dimension,Soil_parameters)
(material_vector,gamma_vector,Permeability_curve_Perm_vector,Shrinkage_curve_vector,Water_retention_curve_vector, Compression_vector,Emission_vector) = soilbuilder(Soil_dimension,Soil_parameters,Soil_dimension.node_vector,Soil_dimension.internode_vector)


S_e_ini = Soil_states.theta/Soil_states.n_ini


Soil_dimension.node_initial=get_node_distance(Soil_dimension)

" ### RUN PROGRAM ###" 

(Soil_stresses) = calc_stress_ini_h (Soil_dimension,gamma_vector,Soil_states,Water_retention_curve_vector,Soil_stresses)

Soil_stresses.sigma_effective_node_ini = np.array(Soil_stresses.sigma_effective_node)
Soil_states.Preconsolidation_pressure =  Soil_stresses.sigma_effective_node_ini*Compression_vector.OCR
#--------------------------------------------------------
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
#------------------------------------------------------------------------
#------------------------------------------------------------------------
#Initialize tracker ------------------------------
sigma_tracker = SigmaTracker()
#------------------------------------------------------------------------
#------------------------------------------------------------------------

#@profile
def equations(simulation_time, in_in, sigma_tracker):
    global gamma_vector
    global Compression_vector
    global Shrinkage_curve_vector
    global Water_retention_curve_vector
    global Permeability_cruve_Perm_vector
    global Emission_vector
    global Soil_dimension
    global Bctop
    global Bcbottom
    
    #print("~~~~~~~~~~~~~~~~~~")
    print ("time", "%.2f" % simulation_time)
    "~~~~~~~setup_local_parameters~~~~~~"    
    
    Shrinkage_curve_node = Shrinkage_curve_vector

    Water_retention_curve_node = Water_retention_curve_vector
    "~~~~ retreive parameters from loop~~~~~~~"
    
    " here the / 4 means that there are three sets of parameters changing in dx/dt"
    
    nr_elements = int(np.size(in_in)/6)  


    Pressure,theta,e = untangle_parameters(in_in.copy(),3)
    Soil_states.theta = theta
    Soil_states.e = e
    n = e/(1+e)
    S_e = np.clip(update_S_e(Water_retention_curve_node, theta), 0.0, 1.0)
    saturated = (S_e >= 1).astype(int)
    unsaturated = (S_e < 1).astype(int)  
      
    "~~~~ Set boundary ~~~~~"
    
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
    
    #print(Bctop_in.pressure)
                                           # kPa
    Bctop_in.devaporation = Bctop.devaporation[t_interval]
    Bctop_in.trans = Bctop.trans[t_interval]

    t_interval = len(Bcbottom.time)-1
    for i in range(len(Bcbottom.time)):
        if simulation_time < Bcbottom.time[i]:
            t_interval = t_interval-1  


    Bcbottom_in = Bc_bottom()
    Bcbottom_in.pressure = Bcbottom.pressure[t_interval]                                                         # kPa
    
    
    "~~~~~~~Update soil_dimension~~~~~~"
    Soil_dimension.node_vector = get_node_coordinates(Soil_dimension)
    Soil_dimension.node_distance = get_node_distance(Soil_dimension)  
    Soil_dimension.node_distance = (1-((Soil_states.e_ini-e)/(1+Soil_states.e_ini)))*Soil_dimension.node_initial
         
    "~~~~~~Update_soil_stresses~~~~~~~"
    sigma_total = np.zeros(len(e))
    sigma_total = update_sigma_total (Soil_dimension.node_distance,gamma_vector,Soil_states,Water_retention_curve_vector,Bctop_in)
    sigma_eff = np.abs(sigma_total - Pressure*S_e)
    #compute_dsigma_eff_dt------------------------------------------------
    dsigma_eff_dt = sigma_tracker.compute_dsigma_eff_dt(sigma_eff, simulation_time) 
    #---------------------------------------------------------------------

    "Matrix potential "
    psi_m = np.where(Pressure <= 0, Pressure, 0)
    "~~~~ Calculate permeability ~~~~~"
    Soil_calc.K                     = update_K_Kmin (theta,e,S_e,Permeability_curve_Perm_vector,Water_retention_curve_node)
    K= Soil_calc.K/Soil_dimension.node_distance
    "set differentials"
    # First set: pressure 
    d2psi_dz2 = d2f_dz2_non_uniform_h(Pressure,Bctop_in.pressure,Bcbottom_in.pressure,Soil_dimension) #Pressure*0+np.sin(simulation_time)
    dpsi_dt = K*d2psi_dz2

    dS_e_dpsi,dtheta_dpsi  = derrivative_van_genuchten_S_e_h (Water_retention_curve_node, theta,e, psi_m)

    # third set: void ratio
    de_dt,dn_dt = deSSc_dt_h(Compression_vector, Soil_states, S_e, sigma_eff, dsigma_eff_dt,dpsi_dt,Water_retention_curve_vector, Soil_states.Preconsolidation_pressure)
    

    # second set: theta
    dtheta_dt = dpsi_dt*dtheta_dpsi  
      
    
    compression = (de_dt < 0).astype(int)    
    extension = (de_dt > 0).astype(int)  
    
    positive_pressure = (dpsi_dt > 0).astype(int)
    negative_pressure = (dpsi_dt < 0).astype(int)
    

    dpsi_dt = update_dpsi_dt(dpsi_dt, K, compression, saturated, dn_dt)
        
    dtheta_dt, dSat_dt = update_dtheta_and_dSat(dtheta_dt,positive_pressure, extension, saturated, dn_dt,n, theta, unsaturated)


    
    for i in range(len(S_e)):
        if S_e[i] >= 1:
            if dtheta_dt[i] > dn_dt[i] :  
                dtheta_dt[i] = dn_dt[i]    

    out_out = np.concatenate([dpsi_dt,dtheta_dt,de_dt],axis=0)
    #------------------------------------------------------------------------

    return [out_out]
#------------------------------------------------------------------------
#------------------------------------------------------------------------
"Set input" 
# First set is pressure #
input_in = Soil_stresses.psi_ini
# Second set is water content #
input_in = np.append(input_in,Soil_states.theta)
# Third set is void ratio #
input_in = np.append(input_in,Soil_states.e_ini) 


"This is the code to get the numerical differential equation"

start = time.perf_counter() # measure the performance
sol = solve_ivp(equations, [simulation_time[0], simulation_time[-1]],input_in,method='RK23',t_eval=simulation_time,max_step = max_time_step, rtol = rtol ,args=(sigma_tracker,))#, atol = atol)

end = time.perf_counter()# measure the performance
print("Elapsed  = {}s".format((end - start)))
print("Success:", sol.success)
print("Message:", sol.message)
print("Final time reached:", sol.t[-1])


toc = timeit.default_timer()  



"Present the results - get the results"

nr_elements = int(np.size(input_in)/3)
    # First set is pressure #
Pressure_results = sol.y[0:nr_elements]
    # Second set is water content #
theta_results = sol.y[nr_elements:2*nr_elements]
    # Third set is void ratio #
e_results = sol.y[2*nr_elements:3*nr_elements]
" Extra "




#print (Soil_states.e_ini)
"calculate settlement"

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
#-------------------------------------------------------------
" ~~~~~~~~~~~ "
n_results = e_results/(1+e_results)
#Se_results = theta_results/n_results
Se_results = theta_results/e_results

#AU_results = 0*Emission_results

#AU_results = ((0.05*(Temperature_results-(-10)))**2)*sci.stats.beta.pdf(Se_results,2.59,1.84)


#----------------
#output
##----------------
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

#--------
" ~~~~~~~~~~~~ "
"Plot results"

plt.figure()
plt.plot(date_list,(colum_heigth+Soil_dimension.Depth_bottom)*1000)
plt.xlabel('time [days]')
plt.ylabel('displacement [mm]')
plt.title('displacement')
plt.xticks(rotation=45)
plt.subplots_adjust(left=0.2, bottom=0.2)
plt.savefig('../results/figure/displacement.png')

#------------
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
