# -*- coding: utf-8 -*-
"""
Created on Fri Jul 23 14:27:28 2021

@author: tomdegast

This model links the change in effective stress due to transient flow to a consitutitve modelincluding creep (Modelled after Soft Soil creep). 
The permeability is changed based on watercontent and void ratio. 
Additionally this model can be used to calculate GHG-production rated based on saturation, temperature and available source material


"""
"The following function libraries are used in the programme and are local"
" ## local files##" 

#from main_program import *
from build_domain import *
from NSDWLs_errors import *
from NSDWLs_states import *
from NSDWLs_equations import *
from NSDWLs_potential import *

#from NSDWLs_calc import *
#from NSDWLs_flux import *

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
parent_dir = os.path.abspath("..")
sys.path.insert(0, parent_dir)
from Input_soil_material_parameters import *
from input import *
sys.path.pop(0)
os.chdir(current_dir)
#------

#input file-----------------------------
infile_tp = "input/top_pressure/top_pressure_swp.dat"
#infile_bp = "input/bottom_pressure/bottom_pressure_swp.dat"
infile_wc_ini = "input/wc_ini/wc_ini_swp.dat"
infile_p_ini = "input/pressure_ini/pressure_ini_swp.dat"
df_tp = pd.read_csv(infile_tp, delim_whitespace=True, skiprows=1, names=["date", "pressure"])
#df_bp = pd.read_csv(infile_bp, delim_whitespace=True, skiprows=1, names=["date", "pressure"])
df_wc_ini = pd.read_csv(infile_wc_ini, delim_whitespace=True, skiprows=1, names=["depth", "wcontent"])
df_pressure_ini = pd.read_csv(infile_p_ini, delim_whitespace=True, skiprows=1, names=["depth", "pressure"])
#---------------------------------------
#simulaiton period----------------------
TSTART = tstart
TEND = tend
            

" #### input parameters "
tic = timeit.default_timer()

" Initiation of he soil domain, give top, bottom and precision"
" precision directly influences calculation time. If the precision is small, the relative error required for accuracy will increase calulation time"
Soil_dimension = Soil_domain()
Soil_dimension.Depth_top = 0.00                                                    # m size
Soil_dimension.Depth_bottom = depths[-1]                                               # m size
Soil_dimension.precision = 0.1                                               # m node size

" Layers are added to the domain."
" Soil properties are collected in the library Input_parameters_definitiions "

Soil_dimension.Depth_o_soils = depths                # bottom of soil layer
#Soil_dimension.soils = ["Test_soil", "Test_soil", "Test_soil"]                                      # soil types used
Soil_dimension.soils = types                                      # soil types used

""
print(Soil_dimension.Depth_bottom)
print(Soil_dimension.Depth_o_soils)
print(Soil_dimension.soils)

Soil_dimension.Number_o_soils = len(Soil_dimension.Depth_o_soils)                                                 # number of soils 
Soil_dimension.domainsize = Soil_dimension.Depth_top-Soil_dimension.Depth_bottom




#TEST MATERIALS BELOW
#Clay
#Peat

#Soil_parameters = Soil_parameters()

#Soil_states = Soil_states()
Soil_states = Soil_states()

#Soil_stresses = Soil_stresses()
Soil_stresses = Soil_stresses()

Soil_calc = Soil_calc()
swp = swp()

" ## Set calculation time, simulation steps to be reported, maximum time step and relative tolerance ##"
start_date = datetime.strptime(TSTART, "%Y-%m-%d")
end_date = datetime.strptime(TEND, "%Y-%m-%d")
minimum_accuracy_time_step = 1               #days  
max_time_step = 1.0  
rtol =1e-4 #relative tolerance
atol =1e-16 #relative tolerance

end_time = (end_date - start_date).days
simssteps = (end_date - start_date).days
date_list = [start_date + timedelta(days=i) for i in range((end_date - start_date).days + 1)]
simulation_time =  np.linspace (0,end_time,simssteps+1)



"## Set boundary conditions"

"Top boundary conditions"
Bctop = Bc_top()
Bctop.time =     [0, 50]
Bctop.pressure = [-13.0375, 0]
#Initial conditoin from SWAP----------
#Bctop.time = df_tp["date"].tolist()
#Bctop.pressure = df_tp["pressure"].tolist()   
#-------------------------------------
#Bctop.trans =    [1,  10,   1,   1,   1,  1,   1, 1 ,  1,   1,   1]
Bctop.trans = np.ones(len(Bctop.pressure))
Bctop.Temp = np.ones(len(Bctop.pressure))*10
Bctop.load_increments = [1.9809, 3.0183, 1.9867, 2.992, 4.9896, 5.0048, 19.8894, -19.8379, 19.9315, 40.0327, 59.9423]
Bctop.load = np.cumsum(Bctop.load_increments)
Bctop.time_incr = [0, 1.00393, 8.02939, 9.03358, 10.0375, 11.04165,12.04558119, 13.05007096, 14.05408645, 15.05812487, 22.08338878]

"Top boundary conditions in equation form <-> the choice to use an equation or steps is made further on in the code"
#Make a boundary condition
def Boundary_equation(time):
    t = time 
    P = 5.0*np.sin(np.pi*t*(1/100)+0)-10.0*np.sin(np.pi*t*(1/50)+1)-2.5*np.sin(np.pi*t*(1/500))
    T = np.sin(t*(6/365)-1.8)+7.5
    #P=-1.25
    return (P,T)

"Bottom boundary conditions"        
                               # kPa
Bctop.devaporation = np.zeros(len(Bctop.pressure))

Bcbottom = Bc_bottom()
Bcbottom.time = [0]
Bcbottom.pressure = [49]
#Bcbottom.time = df_bp["date"].tolist()
#Bcbottom.pressure = df_bp["pressure"].tolist()  
Bcbottom.Temp = np.ones(len(Bcbottom.pressure))*10                                                       # kPa
Bcbottom.dflow = np.zeros(len(Bcbottom.pressure))
Bcbottom.trans = np.ones(len(Bcbottom.pressure))

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
Soil_states.n_ini = Soil_states.e_ini/(1+Soil_states.e_ini)
#Soil_states.theta_ini[0]=Soil_states.theta_ini[0]*0.9
Soil_states.theta = Soil_states.theta_ini
#Soil_states.e_ini = Soil_states.theta
Soil_dimension.internode_vector = modeldomain (Soil_dimension,Soil_parameters)
(material_vector,gamma_vector,Permeability_curve_Perm_vector,Shrinkage_curve_vector,Water_retention_curve_vector, Compression_vector,Emission_vector) = soilbuilder(Soil_dimension,Soil_parameters,Soil_dimension.node_vector,Soil_dimension.internode_vector)

Temperature = np.full_like(Soil_states.theta, 10)

S_e_ini = Soil_states.theta_ini/Soil_states.n_ini

#required for GHG emission
(Emission_ini,dEmission_ini) = respiration_rate_curve (S_e_ini, Soil_states.e, Soil_states.theta, np.zeros_like(Soil_states.theta), 0*Temperature,Temperature, Water_retention_curve_vector,Emission_vector)

Soil_dimension.node_initial=get_node_distance(Soil_dimension)

" ### RUN PROGRAM ###" 

(Soil_stresses) = calc_stress (Soil_dimension,gamma_vector,Soil_states,Water_retention_curve_vector,Soil_stresses)

Soil_stresses.sigma_effective_node_ini = np.array(Soil_stresses.sigma_effective_node)
Soil_states.Preconsolidation_pressure =  Soil_stresses.sigma_effective_node_ini*Compression_vector.OCR

#--------------------------------------------------------
#Define a wrapper function to store previous values------
class SigmaTracker:
    def __init__(self):
        self.prev_sigma_eff = None  # Store previous sigma value
        self.prev_t = None      # Store previous time

    def compute_dsigma_eff_dt(self, sigma_eff, t):
        if self.prev_sigma_eff is None:  # First time step, no derivative yet
            dsigma_eff_dt = 0  # Assume initial derivative is 0
        else:
            dt = t - self.prev_t  # Compute time step
            if dt == 0:
                dsigma_eff_dt = 0  # Avoid division by zero
            else:
                dsigma_eff_dt = (sigma_eff - self.prev_sigma_eff) / dt  # Finite difference

        # Update stored values for next step
        self.prev_sigma_eff = sigma_eff
        self.prev_t = t
        return dsigma_eff_dt
#------------------------------------------------------------------------
#------------------------------------------------------------------------
#Initialize tracker ------------------------------
sigma_tracker = SigmaTracker()
#------------------------------------------------------------------------
#------------------------------------------------------------------------
def equations(simulation_time, in_in, sigma_tracker):
    #global time_step
    #global Soil_calc
    #global Results_states
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
    #print(in_in[0:nr_elements])
    #print(in_in[nr_elements:2*nr_elements])
    #print(in_in[2*nr_elements:3*nr_elements])
    #print(in_in[3*nr_elements:4*nr_elements])
    #print(in_in[4*nr_elements:5*nr_elements])
    #print(in_in[5*nr_elements:6*nr_elements])

    Pressure,theta,e,Preconsolidation_pressure,T,emm = untangle_parameters(in_in.copy(),6)
    Soil_states.theta = theta
    Soil_states.e = e
    n = e/(1+e)
    #S_e                   = theta/n #update_S_e(Water_retention_curve_node,theta)
    #S_e                   = update_S_e(Water_retention_curve_node,theta) 
    S_e = np.clip(update_S_e(Water_retention_curve_node, theta), 0.0, 1.0)
    saturated = (S_e >= 1).astype(int)
    unsaturated = (S_e < 1).astype(int)
    #print(saturated)
    #print("Se","%.1f" % S_e[3])    
      
    "~~~~ Set boundary ~~~~~"
    
    t_interval = len(Bctop.time)-1
    for i in range(len(Bctop.time)):
        #print(i)
        if simulation_time < Bctop.time[i]:
            t_interval = t_interval-1         

    Bctop_in = Bc_top()
    Bctop_in.pressure =  Bctop.pressure[t_interval]
    Bctop_in.Temp = Bctop.Temp[t_interval]

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
    Bcbottom_in.pressure = Bcbottom.pressure[t_interval] 
    Bcbottom_in.Temp = Bcbottom.Temp[t_interval]                                                         # kPa
    Bcbottom_in.dflow = Bcbottom.dflow[t_interval]
    Bcbottom_in.trans = Bcbottom.trans[t_interval]
    
    #Bctop_in.pressure,Bctop_in.Temp = Boundary_equation(simulation_time)
    
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
    #print(Pressure[0])
    #print(sigma_total[0])
    #print(sigma_eff[0])

    "Matrix potential "
    psi_m = np.where(Pressure <= 0, Pressure, 0)
    
    "~~~~ Calculate permeability ~~~~~"
    Soil_calc.K                     = update_K_h (theta,e,S_e,Permeability_curve_Perm_vector,Water_retention_curve_node)
    K= Soil_calc.K/Soil_dimension.node_distance
    
    "set differentials"
    # First set is pressure 
    d2psi_dz2 = d2f_dz2_non_uniform_h(Pressure,Bctop_in.pressure,Bcbottom_in.pressure,Soil_dimension) #Pressure*0+np.sin(simulation_time)
    dpsi_dt = K*d2psi_dz2
    
    #print("Pressure", Pressure[0:4])
    #dS_e_dpsi = derrivative_van_genuchten_S_e (Water_retention_curve_node, theta,e, psi_m)
    #dtheta_dpsi = dS_e_dpsi*n

    dS_e_dpsi,dtheta_dpsi  = derrivative_van_genuchten_S_e_h (Water_retention_curve_node, theta,e, psi_m)

    # third set is void ratio
    de_dt,dn_dt,dPP_dt = deSSc_dt_h(Compression_vector, Soil_states, S_e, sigma_eff, dsigma_eff_dt,dpsi_dt,Water_retention_curve_vector, Preconsolidation_pressure)
    
    
    d2T_dz2 = d2f_dz2_non_uniform_h(T,Bctop_in.Temp,Bcbottom_in.Temp,Soil_dimension)   
    dT_dt = 0.044064*d2T_dz2     

    # second set is theta

    dtheta_dt = dpsi_dt*dtheta_dpsi  
      
    
    compression = (de_dt < 0).astype(int)    
    extension = (de_dt > 0).astype(int)  
    
    positive_pressure = (dpsi_dt > 0).astype(int)
    negative_pressure = (dpsi_dt < 0).astype(int)
    
  
    dpsi_dt = dpsi_dt + (-1/K)*compression*saturated*dn_dt
    dpsi_dt = np.where(np.isfinite(dpsi_dt), dpsi_dt, 0.0)        
  
    #print((positive_pressure*extension-1)*-1)
    #print((positive_pressure*extension))
    
    dtheta_dt = -(dtheta_dt - positive_pressure*extension*saturated*(dn_dt) )
    
    for i in range(len(S_e)):
        if S_e[i] >= 1:
            if dtheta_dt[i] > dn_dt[i] :  
                dtheta_dt[i] = dn_dt[i]    

    #dSat_dt = dtheta_dt/n - unsaturated*dn_dt/n
    dSat_dt = (n*dn_dt-theta*dn_dt)/n**2 - unsaturated*dn_dt/n
    #dSat_dt = dtheta_dt/e - de_dt/e

    # emission
    
    #(Emission,dEmission) = respiration_rate_curve (S_e,e, theta,dSat_dt, dT_dt, T, Water_retention_curve_vector,Emission_vector)

    demm_dt = 0*dSat_dt
    
    for i in range(len(emm)):
        if emm[i] <= 0 : 
            if demm_dt[i] <= 0:
                demm_dt[i] = 0
                
    out_out = np.concatenate([dpsi_dt,dtheta_dt,de_dt,dPP_dt,dT_dt,demm_dt],axis=0)
    #------------------------------------------------------------------------
    if not np.all(np.isfinite(dSat_dt)):
        print("NaN or Inf in dSat_dt!", dSat_dt)
        raise ValueError("Non-finite derivative encountered!")

    if not np.all(np.isfinite(S_e)):
        print("NaN or Inf in S_e!", S_e)
        raise ValueError("Non-finite derivative encountered!")

    if not np.all(np.isfinite(K)):
        print("NaN or Inf in K!", K)
        raise ValueError("Non-finite derivative encountered!")

    if not np.all(np.isfinite(dpsi_dt)):
        print("NaN or Inf in dpsi_dt!", dpsi_dt)
        raise ValueError("Non-finite derivative encountered!")

    if not np.all(np.isfinite(dtheta_dt)):
        print("NaN or Inf in dtheta_dt!", dtheta_dt)
        raise ValueError("Non-finite derivative encountered!")

    if not np.all(np.isfinite(de_dt)):
        print("NaN or Inf in de_dt!", de_dt)
        raise ValueError("Non-finite derivative encountered!")

    if not np.all(np.isfinite(dPP_dt)):
        print("NaN or Inf in dPP_dt!", dPP_dt)
        raise ValueError("Non-finite derivative encountered!")
    #------------------------------------------------------------------------
    return [out_out]
#------------------------------------------------------------------------
#------------------------------------------------------------------------
"Set input" 
# First set is pressure #
input_in = Soil_stresses.psi+Soil_stresses.sigma_pressure_sum_node
# Second set is water content #
input_in = np.append(input_in,Soil_states.theta)
# Third set is void ratio #
input_in = np.append(input_in,Soil_states.e_ini) 
# Fourth set is max effective stress #
input_in = np.append(input_in,Soil_states.Preconsolidation_pressure)
# Fifth set is Temperature#
input_in = np.append(input_in,Temperature)
# Sixth set is emmision#
input_in = np.append(input_in,Emission_ini)

"This is the code to get the numerical differential equation"
sol = solve_ivp(equations, [simulation_time[0], simulation_time[-1]],input_in,method='DOP853',t_eval=simulation_time,max_step = max_time_step, rtol = rtol ,args=(sigma_tracker,))#, atol = atol)

#sol = solve_ivp(equations, [simulation_time[0], simulation_time[-1]],input_in,method='RK45',t_eval=simulation_time,max_step = max_time_step, rtol = rtol)#, atol = atol)
#sol = solve_ivp(equations, [simulation_time[0], simulation_time[-1]],input_in,method='LSODA',t_eval=simulation_time)
#sol = solve_ivp(equations, [simulation_time[0], simulation_time[-1]],input_in,method='DOP853',t_eval=simulation_time)
print("Success:", sol.success)
print("Message:", sol.message)
print("Final time reached:", sol.t[-1])


toc = timeit.default_timer()  



"Present the results - get the results"

nr_elements = int(np.size(input_in)/6)
    # First set is pressure #
Pressure_results = sol.y[0:nr_elements]
    # Second set is water content #
theta_results = sol.y[nr_elements:2*nr_elements]
    # Third set is void ratio #
e_results = sol.y[2*nr_elements:3*nr_elements]
    # Fourth set is max effective stress #
Preconsolidation_pressure_results = sol.y[3*nr_elements:4*nr_elements]
    # Fifth set is Temperature#
Temperature_results = sol.y[4*nr_elements:5*nr_elements]
    # Fifth set is Emmision#
Emission_results = sol.y[5*nr_elements:6*nr_elements]
" Extra "


#diffrence in pressure
dPressure_results = 0*Pressure_results
dtheta_results = 0*theta_results
de_results = 0*e_results


for i in range(len(Pressure_results[0,:])):
    dPressure_results[:,i] = Pressure_results[:,i]-Pressure_results[:,0]
    dtheta_results[:,i] = theta_results[:,i]-theta_results[:,0]
    de_results[:,i] = e_results[:,i]-e_results[:,0]

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

pressure_top = sol.t*0
Temperature_top =  sol.t*0

for i in range(len(sol.t)):
    pressure_top[i], Temperature_top[i] = Boundary_equation(sol.t[i])

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

AU_results = 0*Emission_results

AU_results = ((0.05*(Temperature_results-(-10)))**2)*sci.stats.beta.pdf(Se_results,2.59,1.84)


#----------------
#output
##----------------
with open('output/result/displacement.dat', 'w') as dis:
    dis.write(f"{'':<2}{'time[days]':<15}  {'displacement [mm]':<15}\n")
    for i in range(len(date_list)):
        dis.write(f"{'':<2}{date_list[i].strftime('%Y-%m-%d'):<15}  {(colum_heigth[i]+Soil_dimension.Depth_bottom)*1000:<15}\n")

with open('output/result/Pressure_top.dat', 'w') as dis:
    dis.write(f"{'':<2}{'time[days]':<15} {'Pressure top':<15} \n")
    for i in range(len(date_list)):
        dis.write(f"{'':<2}{date_list[i].strftime('%Y-%m-%d'):<15} {Pressure_results[0,i]:<15} \n")

with open('output/result/Se.dat', 'w') as dis:
    dis.write(f"{'':<2}{'time[days]':<15} {'depth':<15} {'Se':<15}\n")
    for i in range(len(date_list)):
        dis.write(f"{'':<2}{date_list[i].strftime('%Y-%m-%d'):<15} {colum_node_heigth[1,i]+Soil_dimension.Depth_bottom:<15} {theta_results[1,i]*(1+e_results[1,i])/e_results[1,i]:<15}\n")

with open('output/result/theta.dat', 'w') as dis:
    dis.write(f"{'':<2}{'time[days]':<15} {'depth':<15}  {'theta':<15}\n")
    for i in range(len(date_list)):
        dis.write(f"{'':<2}{date_list[i].strftime('%Y-%m-%d'):<15} {colum_node_heigth[1,i]+Soil_dimension.Depth_bottom:<15} {theta_results[1,i]:<15}\n")

with open('output/result/n.dat', 'w') as dis:
    dis.write(f"{'':<2}{'time[days]':<15}  {'depth':<15} {'n':<15}\n")
    for i in range(len(date_list)):
        dis.write(f"{'':<2}{date_list[i].strftime('%Y-%m-%d'):<15} {colum_node_heigth[1,i]+Soil_dimension.Depth_bottom:<15} {e_results[1,i]/(1+e_results[1,i]):<15}\n")


#--------
" ~~~~~~~~~~~~ "
"Plot results"

#plt.figure()
#plt.plot(date_list,Bctop.pressure)
#plt.xlabel('time [days]')
#plt.ylabel('pressure_top [kPa]')
#plt.title('Pressure')
#plt.xticks(rotation=45)
#plt.savefig('output/pressure_top_input.png')
plt.figure()
plt.plot(simulation_time,Pressure_results[0,:])
plt.xlabel('time [days]')
plt.ylabel('pressure_top [kPa]')
plt.title('Pressure')
plt.savefig('output/pressure_top_0.png')

plt.figure()
plt.plot(simulation_time,(colum_heigth+Soil_dimension.Depth_bottom)*1000)
plt.xlabel('time [days]')
plt.ylabel('displacement [mm]')
plt.title('displacement')
plt.xticks(rotation=45)
plt.subplots_adjust(left=0.2, bottom=0.2)
plt.savefig('output/displacement.png')

plt.figure()
for i in range(len(e_results[:,0])):  
    plt.plot(simulation_time,(colum_node_heigth[i,:]+Soil_dimension.Depth_bottom))
plt.xlabel('time [days]')
plt.ylabel('displacement [m]')
plt.title('displacement')
plt.xticks(rotation=45)
plt.subplots_adjust(left=0.2, bottom=0.2)
plt.savefig('output/displacement_layers.png')


#------------
anc =[]
plt.figure()
anchor_to_mv = [-7.80, -5.30, -1.95, -1.30, -0.50, -0.07] # [m]
filtered_anchor_to_mv = [value for value in anchor_to_mv if value >= Soil_dimension.Depth_bottom]
i_anchor = 0
anchor_depth = np.zeros((len(filtered_anchor_to_mv), len(simulation_time)))
for i in range(len(e_results[:,0])):
    if i_anchor < len(filtered_anchor_to_mv):
        if colum_node_heigth[i,0]+Soil_dimension.Depth_bottom <= filtered_anchor_to_mv[-(i_anchor+1)]:
            anchor_depth[i_anchor,:]=colum_node_heigth[i,:]
            i_anchor = i_anchor+1
            plt.plot(simulation_time,colum_node_heigth[i,:]+Soil_dimension.Depth_bottom)
            anc.append(i)

plt.xlabel('time [days]')
plt.ylabel('displacement [m]')
plt.title('displacement at the anchor depth')
plt.xticks(rotation=45)
plt.subplots_adjust(left=0.2, bottom=0.2)
plt.savefig('output/displacement_anchor.png')
#------------
#------------
plt.figure()
plt.plot(simulation_time, (colum_node_heigth[anc[1], :] + Soil_dimension.Depth_bottom) * 1000)
plt.xlabel('Time [days]')
plt.ylabel('Displacement [mm]')
plt.title('Displacement at anchor 2')
plt.xticks(rotation=45)
ax = plt.gca()
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
plt.tight_layout()
plt.savefig('output/displacement_anchor_2.png')
plt.close() 


plt.figure()
plt.plot(simulation_time, (colum_node_heigth[anc[2], :] + Soil_dimension.Depth_bottom) * 1000)
plt.xlabel('Time [days]')
plt.ylabel('Displacement [mm]')
plt.title('Displacement at anchor 3')
plt.xticks(rotation=45)
ax = plt.gca()
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
plt.tight_layout()
plt.savefig('output/displacement_anchor_3.png')
plt.close()


plt.figure()
plt.plot(simulation_time, (colum_node_heigth[anc[3], :] + Soil_dimension.Depth_bottom) * 1000)
plt.xlabel('Time [days]')
plt.ylabel('Displacement [mm]')
plt.title('Displacement at anchor 4')
plt.xticks(rotation=45)
ax = plt.gca()
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
plt.tight_layout()
plt.savefig('output/displacement_anchor_4.png')
plt.close()
#------------
#------------
#plt.figure()
#i_str = 1
#strain = []
#legend_labels = []
#org_h = []
#for i in range(len(filtered_anchor_to_mv)-1):
#    org_h.append(filtered_anchor_to_mv[i+1] - filtered_anchor_to_mv[i])
#
#for i in range(len(anchor_depth)-1):
#    plt.plot(date_list,((anchor_depth[i,:] - anchor_depth[i+1,:])-org_h[len(org_h)-1-i])/org_h[len(org_h)-1-i])
#    legend_labels.append(f'Zone {len(anchor_to_mv)-len(filtered_anchor_to_mv)+2+i}')
#
#plt.legend(legend_labels, loc='upper right') 
#plt.xlabel('time [days]')
#plt.ylabel('strain [m/m]')
#plt.title('Strain between anchors')
#plt.xticks(rotation=45)
#plt.subplots_adjust(left=0.2, bottom=0.2)  
#plt.savefig('output/displacement_strain.png')

#------------
#------------
abs_max = np.max(np.abs(Pressure_results))

plt.figure()
plt.pcolormesh(simulation_time, Soil_dimension.node_vector, Pressure_results,
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
plt.savefig('output/pressure.png')
#--------------

#plt.figure()
#plt.pcolormesh(sol.t,Soil_dimension.node_vector,dPressure_results, cmap='RdBu', vmin=np.min(dPressure_results), vmax=np.max(dPressure_results))
##plt.pcolormesh(sol.t,Soil_dimension.node_vector,Pressure_results, cmap='RdBu', vmin=-10, vmax=100)
#for i in range(len(Soil_dimension.Depth_o_soils)-1):
#    plt.axhline(y=Soil_dimension.Depth_o_soils[i], color='k', linestyle='--')
#plt.title("delta pressure")
#plt.xlabel("time [days]")
#plt.ylabel("depth [m]")
#plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
#cax = plt.axes((0.85,0.1,0.075,0.8))
#plt.colorbar(cax=cax)
#plt.title("kPa")
#plt.savefig('output/delta_pressure.png')


#plt.figure()
#plt.pcolormesh(sol.t,Soil_dimension.node_vector,de_results, cmap='RdBu', vmin=np.min(de_results), vmax=np.max(de_results))
#plt.pcolormesh(sol.t,Soil_dimension.node_vector,Pressure_results, cmap='RdBu', vmin=-10, vmax=100)
#for i in range(len(Soil_dimension.Depth_o_soils)-1):
#    plt.axhline(y=Soil_dimension.Depth_o_soils[i], color='k', linestyle='--')
#plt.title("delta e")
#plt.xlabel("time [days]")
#plt.ylabel("depth [m]")
#plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
#cax = plt.axes((0.85,0.1,0.075,0.8))
#plt.colorbar(cax=cax)
#plt.title("-")
#plt.savefig('output/delta_e.png')

#plt.figure()
#norm = mcolors.Normalize(vmin=np.min(dtheta_results-de_results), vmax=np.max(dtheta_results-de_results))
#plt.pcolormesh(sol.t,Soil_dimension.node_vector,dtheta_results-de_results, cmap='Blue', norm=norm)
##plt.pcolormesh(sol.t,Soil_dimension.node_vector,Pressure_results, cmap='RdBu', vmin=-10, vmax=100)
#for i in range(len(Soil_dimension.Depth_o_soils)-1):
#    plt.axhline(y=Soil_dimension.Depth_o_soils[i], color='k', linestyle='--')
#plt.title("delta theta-de")
#plt.xlabel("time [days]")
#plt.ylabel("depth [m]")
#plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
#cax = plt.axes((0.85,0.1,0.075,0.8))
#plt.colorbar(cax=cax)
#plt.title("-")
#plt.savefig('output/delta_theta_de.png')

plt.figure()
#plt.pcolormesh(sol.t,Soil_dimension.node_vector,theta_results, cmap='RdBu', vmin=np.min(theta_results), vmax=np.max(theta_results))
plt.pcolormesh(simulation_time,Soil_dimension.node_vector,theta_results, cmap='Blues', vmin=np.min(theta_results), vmax=np.max(theta_results))
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
plt.savefig('output/theta.png')


plt.figure()
vG_WCR = Water_retention_curve_vector.vG_WCR[:, np.newaxis]  
vG_WCS = Water_retention_curve_vector.vG_WCS[:, np.newaxis]
norm = mcolors.Normalize(vmin=np.min((theta_results-vG_WCR)/(vG_WCS-vG_WCR)), vmax=np.max((theta_results-vG_WCR)/(vG_WCS-vG_WCR)))
plt.pcolormesh(simulation_time,Soil_dimension.node_vector,(theta_results-vG_WCR)/(vG_WCS-vG_WCR), cmap='Blues', norm=norm)
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
plt.savefig('output/Se_vG.png')


plt.figure()
norm = mcolors.Normalize(vmin=np.min(e_results/(1+e_results)), vmax=np.max(e_results/(1+e_results)))
plt.pcolormesh(simulation_time, Soil_dimension.node_vector, e_results/(1+e_results), cmap='Blues', norm=norm)
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
plt.savefig('output/porosity.png')

#plt.figure()
#plt.pcolormesh(date_list,Soil_dimension.node_vector,Preconsolidation_pressure_results, cmap='cividis')
#for i in range(len(Soil_dimension.Depth_o_soils)-1):
#    plt.axhline(y=Soil_dimension.Depth_o_soils[i], color='k', linestyle='--')
#plt.title("Plastic point")
#plt.xlabel("time [days]")
#plt.ylabel("depth [m]")
#plt.xticks(rotation=45)
#plt.subplots_adjust(left=0.2, bottom=0.2) 
#plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
#cax = plt.axes((0.85,0.1,0.075,0.8))
#plt.colorbar(cax=cax)
#plt.title("[-]")
#plt.savefig('output/Preconsolidation_pressure.png')


#plt.figure()
#plt.pcolormesh(date_list,Soil_dimension.node_vector,Temperature_results, cmap='RdBu')
#for i in range(len(Soil_dimension.Depth_o_soils)-1):
#    plt.axhline(y=Soil_dimension.Depth_o_soils[i], color='k', linestyle='--')
#plt.title("Temperature []")
#plt.xlabel("time [days]")
#plt.ylabel("depth [m]")
#plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
#cax = plt.axes((0.85,0.1,0.075,0.8))
#plt.colorbar(cax=cax)
#plt.title("[-]")
#plt.savefig('output/temperature.png')

#plt.figure()
#plt.pcolormesh(sol.t,Soil_dimension.node_vector,Emission_results, cmap='RdBu')
#for i in range(len(Soil_dimension.Depth_o_soils)-1):
#    plt.axhline(y=Soil_dimension.Depth_o_soils[i], color='k', linestyle='--')
#plt.title("Emission")
#plt.xlabel("time [days]")
#plt.ylabel("depth [m]")
#plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
#cax = plt.axes((0.85,0.1,0.075,0.8))
#plt.colorbar(cax=cax)
#plt.title("[-]")
#plt.savefig('output/emission.png')

#plt.figure()
#plt.pcolormesh(date_list,Soil_dimension.node_vector,AU_results, cmap='RdBu')
#for i in range(len(Soil_dimension.Depth_o_soils)-1):
#    plt.axhline(y=Soil_dimension.Depth_o_soils[i], color='k', linestyle='--')
#plt.title("Emission")
#plt.xlabel("time [days]")
#plt.ylabel("depth [m]")
#plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
#cax = plt.axes((0.85,0.1,0.075,0.8))
#plt.colorbar(cax=cax)
#plt.title("[-]")
#plt.savefig('output/emission_au.png')

#plt.figure()
#plt.pcolormesh(sol.t,Soil_dimension.node_vector,gwl, cmap='RdBu', vmin=0, vmax=5)
#for i in range(len(Soil_dimension.Depth_o_soils)-1):
#    plt.axhline(y=Soil_dimension.Depth_o_soils[i], color='k', linestyle='--')
#plt.title("ground water level")
#plt.xlabel("time [days]")
#plt.ylabel("depth [m]")
#plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
#cax = plt.axes((0.85,0.1,0.075,0.8))
#plt.colorbar(cax=cax)
#plt.title("[-]")
#plt.savefig('output/ground_water_level.png')

plt.figure(figsize=(12, 6))
plt.plot(simulation_time,gwl)
plt.title("Groundwater level")
plt.xlabel("time [days]")
plt.ylabel("depth [m]")
plt.xticks(rotation=45)
plt.subplots_adjust(left=0.2, bottom=0.2) 
plt.savefig('output/ground_water_level.png')

#plt.figure()
#plt.plot(Se_results[:,[0,10,20,30,40,50,60,70,80,90,100]],Soil_dimension.node_vector)
#for i in range(len(Soil_dimension.Depth_o_soils)-1):
#    plt.axhline(y=Soil_dimension.Depth_o_soils[i], color='k', linestyle='--')
##plt.pcolormesh(sol.t,Soil_dimension.node_vector,theta_results, cmap='RdBu', vmin=0, vmax=10)
#plt.title("theta")
#plt.xlabel("Se [-]")
#plt.ylabel("depth [m]")
#plt.savefig('output/theta_Se.png')





#plt.style.use('_mpl-gallery')
#fig, ax = plt.subplots()
#for i in range(Results_states.theta[0,:].size):
#    ax.plot(Results_states.time, Results_states.theta[:,i]-i, linewidth=1.0)
#   
#plt.show()
#
#plt.style.use('_mpl-gallery')
#fig, ax = plt.subplots()
#for i in range(Results_states.theta[0,:].size):
#    ax.plot(Results_states.time, Results_states.e[:,i]-i, linewidth=1.0)
#   
#plt.show()

#sol  = diffusion(sigma_sumsoil_ini, simulation_time, e_ini, Soil_dimension.precision, Soil_dimension.domainsize, Bctop, Bcbottom)

