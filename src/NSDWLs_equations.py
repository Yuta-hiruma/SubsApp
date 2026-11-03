# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 11:36:27 2022

@author: tomdegast
"""
import numpy as np
import scipy as sci
from numba import jit



@jit(forceobj=True, looplift=True,fastmath = True)
def calc_stress_ini_h (Soil_dimension,gamma_vector,Soil_states,Water_retention_curve_vector,Soil_stresses):

    node_size_initial = get_internode_distance(Soil_dimension)
    
    # gamma_vector
    gamma_s     = np.array(gamma_vector.s)
    gamma_w     = np.array(gamma_vector.w)
    
    # WRC
    WCR = np.array(Water_retention_curve_vector.vG_WCR)
    WCS = np.array(Water_retention_curve_vector.vG_WCS)
    
    # States
    theta = np.array(Soil_states.theta)
    e = np.array(Soil_states.e)
    
    # Stresses
    psi_in = Soil_stresses.psi_ini
    # Effective saturation
    S_sat = (theta -WCR) / (WCS -WCR)

    # calculate stotal stress
    sigma_total_node = np.zeros(np.shape(gamma_s))
    sigma_total_node = 9.81*(gamma_s*(1/(e+1))+S_sat*gamma_w*(1-(1/(e+1))))
    
    "calculate sum of total stress from nodes and project on nodes and internodes"
    sigma_total_sum_internode = np.zeros(np.shape(gamma_s))
    sigma_total_sum_node = np.zeros(np.shape(gamma_s))
    sigma_total_sum_internode = np.append(sigma_total_sum_internode,0)

    for i in range(len(sigma_total_node)):
       sigma_total_sum_internode[i+1] = node_size_initial[i]*sigma_total_node[i]+sigma_total_sum_internode[i]    

    for i in range(len(sigma_total_sum_node)):
        sigma_total_sum_node[i] = 0.5 * sigma_total_sum_internode[i] + 0.5 * sigma_total_sum_internode[i+1]


    # step 4 Initial pore pressure
    Soil_stresses.psi = Soil_stresses.psi_ini
    psi = np.array(Soil_stresses.psi)

    # step 5 calculate effective stress
    sigma_effective_node = sigma_total_sum_node - S_sat*psi

    Soil_stresses.sigma_total_node              = np.array(sigma_total_node)
    Soil_stresses.sigma_total_sum_node          = np.array(sigma_total_sum_node)
    Soil_stresses.sigma_total_sum_internode     = np.array(sigma_total_sum_internode)
    Soil_stresses.sigma_effective_node          = np.array(sigma_effective_node)

    return (Soil_stresses)


@jit(forceobj=True, looplift=True,fastmath = True)
def get_node_coordinates(Soil_dimension): 
    internode_vector = Soil_dimension.internode_vector
    node_vector = np.zeros((len(internode_vector)-1))
    for i in range(len(node_vector)):
        node_vector[i] = 0.5*(internode_vector[i]+internode_vector[i+1])
    return (node_vector)

@jit(forceobj=True, looplift=True,fastmath = True)
def get_node_distance(Soil_dimension):
    internode_vector = Soil_dimension.internode_vector
    node_distance = np.zeros(len(internode_vector)-1)
    for i in range(len(node_distance)):
        node_distance[i] = np.abs(internode_vector[i]-internode_vector[i+1])
    return (node_distance)

def get_internode_distance(Soil_dimension):
    internode_vector = Soil_dimension.internode_vector
    internode_distance = np.zeros(len(internode_vector)-1)
    for i in range(len(internode_distance)):
        internode_distance[i] = np.abs(internode_vector[i]-internode_vector[i+1])
    return(internode_distance)

def untangle_parameters (in_in,nr_sets): 
    size_holder = np.size(np.shape(in_in)) 
    if  size_holder> 1:
        in_input = in_in[0]
    else :
        in_input = in_in
         
    nr_elements = int(np.size(in_in)/nr_sets)
    # First set is pressure #
    Pressure = in_input[0:nr_elements]
    # Second set is water content #
    theta = in_input[nr_elements:2*nr_elements]
    # Third set is void ratio #
    e = in_input[2*nr_elements:3*nr_elements]

    return(Pressure, theta, e)



def update_sigma_total (node_size,gamma_vector,Soil_states,Water_retention_curve_vector):
    node_size = node_size
    
    gamma_s     = np.array(gamma_vector.s)
    gamma_w     = np.array(gamma_vector.w)
    
    WCR = np.array(Water_retention_curve_vector.vG_WCR)
    WCS = np.array(Water_retention_curve_vector.vG_WCS)
    
    # States
    theta = np.array(Soil_states.theta)
    e = np.array(Soil_states.e)
    n = e/(1+e)

    
    S_sat = np.clip((theta - WCR) / (WCS - WCR), 0, 1)
    sigma_total_node = np.zeros(np.shape(gamma_s))
    sigma_total_node = 9.81*(gamma_s*(1/(e+1))+S_sat*gamma_w*(1-(1/(e+1))))
    
    sigma_total = np.zeros(len(sigma_total_node))
    sigma_total_in = np.zeros(len(sigma_total_node)+1)
     
    for i in range(len(sigma_total_node)):
        sigma_total_in[i+1] = sigma_total_node[i]*node_size[i]+sigma_total_in[i]
        sigma_total[i] = 0.5*sigma_total_in[i]+0.5*sigma_total_in[i+1]

    return (sigma_total)


def update_K_Kmin (theta,e,S_e,Permeability_curve_Perm_vector,WRC): 
    n = e/(1+e)
    tau = Permeability_curve_Perm_vector.tau
    
    for i in range(len(theta)):
        if theta[i] <= WRC.vG_WCR[i]:
            theta[i] = WRC.vG_WCR[i]
        elif theta[i] >= WRC.vG_WCS[i]:
            theta[i] = WRC.vG_WCS[i]
   
    K_sat = Permeability_curve_Perm_vector.K_tot
    K_min = Permeability_curve_Perm_vector.K_min
    vG_n = WRC.vG_n
    vG_m = WRC.vG_m

    K_out = K_sat*(S_e**tau)*(1-(1-S_e**(1/vG_m)**vG_m))**2
    K_out = np.maximum(K_out, K_min)
    K_out = np.where(S_e < 1,K_out,K_sat)

    return (K_out)

def update_S_e (Water_retention_curve_vector,theta): 
    vG_WCR = Water_retention_curve_vector.vG_WCR
    vG_WCS = Water_retention_curve_vector.vG_WCS
    theta = theta
    S_e = (theta - vG_WCR) / (vG_WCS - vG_WCR)
    return S_e


@jit(forceobj=True, looplift=True,fastmath = True)
def derrivative_van_genuchten_S_e(WRC, theta,e, psi_m):
    Se = np.zeros(np.shape(theta))
    dSe_dpsi = np.zeros(np.shape(theta))
    dtheta_dpsi = np.zeros(np.shape(theta))
    psi_m = -psi_m
    #shrinking
    a = WRC.vG_a
    alpha = WRC.vG_alpha
    n = WRC.vG_n
    m = WRC.vG_m
    theta_s = WRC.vG_WCS
    theta_r = WRC.vG_WCR
    
    Se = (1 + (alpha * np.abs(psi_m))**n)**(-m)
    dSe_dpsi = -m * (1 + (alpha * np.abs(psi_m))**n)**-(m+1) * n * alpha**n * np.abs(psi_m)**(n-1) * np.sign(psi_m)
    dtheta_dpsi = (theta_s - theta_r) * dSe_dpsi

    for i in range(len(psi_m)):
        if psi_m[i] <= 0:
            dSe_dpsi[i] = 0.0

    return dSe_dpsi, dtheta_dpsi



@jit(forceobj=True, looplift=True,fastmath = True)
def d2f_dz2_non_uniform(Parameter,Bctop,Bcbottom,Soil_dimension): 
    """
    2nd derivertive with respect to depth
    """
    Parameter_loc = Parameter
    
    Parameter_loc[0] = Bctop
    Parameter_loc[-1] = Bcbottom
    node_size_loc = 0*Parameter_loc
    d2f_dz2 = 0*Parameter_loc

    
    node_size_loc = Soil_dimension.node_distance
    
    d2f_dz2[0] = (Parameter_loc[2] - 2 * Parameter_loc[1] + Parameter_loc[0]) / (node_size_loc[0]**2)
    d2f_dz2[-1] = (Parameter_loc[-1] - 2 * Parameter_loc[-2] + Parameter_loc[-3]) / (node_size_loc[-2]**2)

    for i in range(1, len(d2f_dz2) - 1):
        dz_forward = node_size_loc[i]
        dz_backward = node_size_loc[i - 1]
        forward_diff = (Parameter_loc[i + 1] - Parameter_loc[i]) / dz_forward
        backward_diff = (Parameter_loc[i] - Parameter_loc[i - 1]) / dz_backward
        d2f_dz2[i] = (forward_diff - backward_diff) / (0.5 * (dz_forward + dz_backward))
    return (d2f_dz2)

def harmonic_mean(a,b):
    c = (2*a*b)/(a+b)
    return c

def logarithmic_mean(a,b):
    if np.isclose(a, b):
        c = a
    else:
        c = (b - a) / np.log(abs(b / a))
    return c


def diffusion(Pressure,Bctop_in,Bcbottom_in,dpsi_top_fun,dpsi_bottom_fun,bc,tp_bc,Soil_dimension,K,t,S_e):
    psi_loc = Pressure
    node_size_loc = 0*psi_loc
    dpsi_dt = 0*psi_loc
    node_size_loc = Soil_dimension.node_distance
    #top boundary condition-------------
    if tp_bc ==0: #Dirichlet 
        psi_loc[0] = Bctop_in.pressure
        if bc == 0:
            dpsi_dt[0] = K[0]*(psi_loc[2] - 2 * psi_loc[1] + psi_loc[0]) / (node_size_loc[0]**2)
        else:
            dpsi_dt[0] = dpsi_top_fun(t)
    else: #Neumann
        k_top = logarithmic_mean(K[0],K[1])
        psi_tp = (Bctop_in.flux*node_size_loc[0]) /k_top + psi_loc[1]
        dpsi_dt[0] =(K[0]*(psi_tp-psi_loc[1]))/node_size_loc[0]

    #bottom boundary condition-------------
    psi_loc[-1] = Bcbottom_in.pressure
    if bc == 0:
        dpsi_dt[-1] = K[-1]*(psi_loc[-1] - 2 * psi_loc[-2] + psi_loc[-3]) / (node_size_loc[-2]**2)
        #xsdpsi_dt[-1] = K[-1]*(psi_loc[-1] -  psi_loc[-2] ) / (node_size_loc[-21]**2)
    else:
        dpsi_dt[-1] = dpsi_bottom_fun(t)        

    for i in range(1, len(dpsi_dt) - 1):
        dz_forward = node_size_loc[i]
        k_forward = logarithmic_mean(K[i-1],K[i])
        dz_backward = node_size_loc[i - 1]
        k_backward = logarithmic_mean(K[i],K[i+1])
        dz = (dz_forward+dz_backward)/2
        dpsi_dt[i] = (k_forward*psi_loc[i - 1]-(k_forward+k_backward)*psi_loc[i]+k_backward*psi_loc[i+1]) / (dz**2)
    return (dpsi_dt)


def deSSc_dt(Compression_vector, Soil_states, S_e, sigma_eff,sigma_eff_p,dsigma_eff_dt, dpsi_dt,Water_retention_curve_vector):
    """
    Soft Soil Creep Model
    """
    e = Soil_states.e
    
    Lambda = Compression_vector.Lambda
    kappa = Compression_vector.Kappa
    mu = Compression_vector.Mu
    Lambda_m = Lambda
    kappa_m = kappa
    mu_m = mu


    Tau_ref = Compression_vector.Tau_ref
    OCR = Compression_vector.OCR

    #dsigma_eff_dt = -dpsi_dt*S_e
    ccott = Compression_vector.ccott
    compressibility = Lambda_m- kappa_m

    exponent = compressibility / mu_m
    #exponent = np.maximum(exponent, 5.0)
    #dPlastic = mu_m * (1 / Tau_ref) * (1/OCR) ** (exponent)
    dPlastic = mu_m * (1 / Tau_ref) * (sigma_eff/sigma_eff_p) ** (exponent)
    dElastic = dsigma_eff_dt*kappa_m/(ccott+sigma_eff)
    deps_dt = -(dElastic + dPlastic) 
    dn_dt = deps_dt/(1+e)**2

    return deps_dt, dn_dt


@jit(forceobj=True, looplift=True,fastmath = True)
def update_dpsi_dt(dpsi_dt, K, compression, saturated, dn_dt):
    dpsi_new = dpsi_dt + (-1 / K) * compression * saturated * dn_dt
    dpsi_new = np.where(np.isfinite(dpsi_new), dpsi_new, 0.0)
    return dpsi_new

@jit(forceobj=True, looplift=True,fastmath = True)
def update_dtheta_and_dSat(dtheta_dt, positive_pressure, extension, saturated, dn_dt,n, theta, unsaturated):#c
    # update dtheta_dt
    dtheta_new = positive_pressure * extension * saturated * dn_dt - dtheta_dt
    # update Sat_dt
    dSat_new = (n * dn_dt - theta * dn_dt) / n**2 - unsaturated * dn_dt / n
    return dtheta_new, dSat_new


def saturated_theta(S_e,dtheta_dt,dn_dt,i):
    if S_e[i] >= 1:
        if dtheta_dt[i] > dn_dt[i] :  
            dtheta_dt[i] = dn_dt[i] 
    return dtheta_dt[i]

def gather(*args):
    return list(args)

#~~~~~~~~~~~~~~classes used~~~~~~~~~~~~~~#
class gamma_vectors:
    def __init__ (self):  
        self.s = []
        self.w = []
        self.o = []
        self.o_ratio = []
        
class Permeability_cruve_Perm_vectors:
    def __init__ (self):  
        self.K_tot = []
        self.K_min = []
        self.tau = []
              
                
class Water_retention_curve_vectors:
    def __init__ (self):  
        self.vG_WCR = []
        self.vG_WCS = []
        self.vG_a = []
        self.vG_alpha = []
        self.vG_m = []
        self.vG_mod = []
        self.vG_n = []


#~~~~~~~~~~~~Results to time writing~~~~~~~~#

class Results_states:
    def __init__(self):
        self.time = []
        self.theta = []
        self.e = []
        self.Se = []
        self.Flow_potential = []
        
class Results_gamma_vectors:
    def __init__ (self):  
        self.time = []
        self.s = []
        self.w = []
        self.o = []
        self.o_ratio = []
        
class Results_Permeability_cruve_Perm_vectors:
    def __init__ (self):  
        self.K_tot = []
        self.K_min = []
        self.tau = []
                
class Results_Water_retention_curve_vectors:
    def __init__ (self):  
        self.time = []
        self.vG_WCR = []
        self.vG_WCS = []
        self.vG_a = []
        self.vG_alpha = []
        self.vG_m = []
        self.vG_mod = []
        self.vG_n = []
        
class Results_Soil_dimentions:
    def __init__ (self):
        self.time = []
        self.internode_distance = []        
