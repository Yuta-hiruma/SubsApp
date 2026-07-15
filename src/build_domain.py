# -*- coding: utf-8 -*-

import numpy as np
import math
from NSDWLs_equations import *

def modeldomain(Soil_dimension,Soil_parameters):
    """
        set model domain (discretizing soil column)
    """
    Depth_top       = Soil_dimension.Depth_top
    Depth_bottom    = Soil_dimension.Depth_bottom
    precision       = Soil_dimension.precision
    Depth_o_soils   = Soil_dimension.Depth_o_soils
    Number_o_soils  = Soil_dimension.Number_o_soils
     
    print ("Top is at", Depth_top, " m", " & Bottom is at", Depth_bottom, " m ")
    
    domain_size = Depth_top-Depth_bottom

    internode_vector = []
    
    internode_vector =np.linspace(Depth_top, Depth_o_soils[0], num=int(np.ceil((Depth_top - Depth_o_soils[0]) /precision)+1))

    if Number_o_soils>1:
        for i in range(len(Depth_o_soils)-1):
            internode_vector = internode_vector[:-1]
            internode_vector = np.append(internode_vector,np.linspace(Depth_o_soils[i], Depth_o_soils[i+1], num=int(np.ceil((Depth_o_soils[i] - Depth_o_soils[i+1]) /precision)+1)))

    print ("the domain is,", domain_size, " m divided in", len(internode_vector)-1, " Nodes.")
    return (internode_vector)


def soilbuilder(Soil_domain,Soil_parameters,node_vector,internode_vector):
    """
        set material properties
    """
    material_vector = []
    i=0
    j=0
    for node in node_vector:
        if node_vector[i] >= Soil_domain.Depth_o_soils[j]:
           material_vector.append (Soil_domain.soils[j])
           i=i+1
        else:
           j=j+1
           material_vector.append (Soil_domain.soils[j]) 
           i=i+1
                
    i=0
    gamma_vector = gamma_vectors()
    gamma_vector.s = np.zeros(np.shape(node_vector))
    gamma_vector.w = np.zeros(np.shape(node_vector))
    
    Permeability_curve_Perm_vector.K_tot = np.zeros(np.shape(node_vector))
    Permeability_curve_Perm_vector.K_min = np.zeros(np.shape(node_vector))
    Permeability_curve_Perm_vector.tau = np.zeros(np.shape(node_vector))
    
    Compression_vector.Kappa = np.zeros(np.shape(node_vector))
    Compression_vector.Lambda = np.zeros(np.shape(node_vector))
    Compression_vector.Mu = np.zeros(np.shape(node_vector))
    Compression_vector.Tau_ref = np.zeros(np.shape(node_vector))
    Compression_vector.OCR = np.zeros(np.shape(node_vector))
    Compression_vector.c = np.zeros(np.shape(node_vector))
    Compression_vector.phi = np.zeros(np.shape(node_vector))
    Compression_vector.ccott = np.zeros(np.shape(node_vector))
    
    Water_retention_curve_vector = Water_retention_curve_vectors()
    Water_retention_curve_vector.vG_WCR = np.zeros(np.shape(node_vector))
    Water_retention_curve_vector.vG_WCS = np.zeros(np.shape(node_vector))
    Water_retention_curve_vector.vG_a = np.zeros(np.shape(node_vector))
    Water_retention_curve_vector.vG_alpha = np.zeros(np.shape(node_vector))
    Water_retention_curve_vector.vG_m = np.zeros(np.shape(node_vector))
    Water_retention_curve_vector.vG_mod = np.zeros(np.shape(node_vector))
    Water_retention_curve_vector.vG_n = np.zeros(np.shape(node_vector))


        
    for material in material_vector:
        gamma_vector.s[i] = Soil_parameters[material]["Material_parameters"]["gamma_s"]
        gamma_vector.w[i] = Soil_parameters[material]["Material_parameters"]["gamma_w"]


        Permeability_curve_Perm_vector.K_tot[i] = Soil_parameters[material]["Permeability_curve_Perm"]["K_tot"]
        Permeability_curve_Perm_vector.K_min[i] = Soil_parameters[material]["Permeability_curve_Perm"]["K_min"]
        Permeability_curve_Perm_vector.tau[i] = Soil_parameters[material]["Permeability_curve_Perm"]["tau"]
        
        Compression_vector.Kappa[i] =Soil_parameters[material]["Compression"]["kappa"]
        Compression_vector.Lambda[i] = Soil_parameters[material]["Compression"]["lambda"]
        Compression_vector.Mu[i] = Soil_parameters[material]["Compression"]["mu"]
        Compression_vector.Tau_ref[i] = Soil_parameters[material]["Compression"]["Tau_ref"]

        Compression_vector.OCR[i] = Soil_parameters[material]["Compression"]["OCR"]
        Compression_vector.c[i] = Soil_parameters[material]["Compression"]["c"]
        Compression_vector.phi[i] = Soil_parameters[material]["Compression"]["phi"]
        Compression_vector.ccott[i] = Compression_vector.c[i]/math.tan(math.radians(Compression_vector.phi[i]))

         
        Water_retention_curve_vector.vG_WCR[i] = Soil_parameters[material]["Water_retention_curve"]["vG_WCR"]
        Water_retention_curve_vector.vG_WCS[i] = Soil_parameters[material]["Water_retention_curve"]["vG_WCS"]
        Water_retention_curve_vector.vG_a[i] = Soil_parameters[material]["Water_retention_curve"]["vG_a"]
        Water_retention_curve_vector.vG_alpha[i] = Soil_parameters[material]["Water_retention_curve"]["vG_alpha"]
        Water_retention_curve_vector.vG_m[i] = Soil_parameters[material]["Water_retention_curve"]["vG_m"]
        Water_retention_curve_vector.vG_mod[i] = Soil_parameters[material]["Water_retention_curve"]["vG_mod"]
        Water_retention_curve_vector.vG_n[i] = Soil_parameters[material]["Water_retention_curve"]["vG_n"]
        
       
        i=i+1
    return (material_vector,gamma_vector,Permeability_curve_Perm_vector,Water_retention_curve_vector,Compression_vector);


def set_ini(Soil_dimension,Soil_parameters,Soil_states,Bctop,swp):
    """
        set initial conditions
    """
    Soil_dimension.internode_vector =  modeldomain(Soil_dimension, Soil_parameters)
    internode_vector = Soil_dimension.internode_vector
    
    Soil_dimension.node_vector = get_node_coordinates(Soil_dimension)
    node_vector = Soil_dimension.node_vector
    
    (material_vector,gamma_vector,Permeability_curve_Perm_vector,Water_retention_curve_vector,Compression_vector) = soilbuilder(Soil_dimension, Soil_parameters,node_vector,internode_vector)
    #Initial water content and pressure from SWAP----------------
    theta_ini = []
    psi_ini = []
    for i in range(len(node_vector)):
        closest_index = min(range(len(swp.depth)), key=lambda j: abs(swp.depth[j] - node_vector[i]))
        theta_ini.append(float(swp.wcontent[closest_index]))  
        if float(swp.wcontent[closest_index]) < 1:
            psi_ini.append(float(swp.pressure[closest_index]))  
        else:
            psi_ini.append(0)
    #----------------------------------------------       
    S_e_ini = np.clip(update_S_e(Water_retention_curve_vector, theta_ini), 0.0, 1.0)
    Soil_states.n = theta_ini/S_e_ini
    Soil_states.e = Soil_states.n /(1-Soil_states.n)
    return (Soil_states,theta_ini,psi_ini)


def effective_stress_ini(gamma_vector, Soil_dimension, Soil_states): 
    """ 
        calculate effective stress based on e_ini 
    """
    gamma_s = gamma_vector.s
    gamma_w = gamma_vector.w
    e = Soil_states.e
    
    sigma_soil = np.zeros(np.shape(gamma_s))
    sigma_soil = 9.81**( gamma_s*(1/(e+1)) + gamma_w*(1-(1/(e+1))))
    sigma_sumsoil = np.zeros(np.shape(gamma_s))
    sigma_sumsoil = np.append(sigma_sumsoil,0)
    node_size = get_internode_distance(Soil_dimension)
    for i in range(len(sigma_soil)):
       sigma_sumsoil[i+1] = node_size[i]*sigma_soil[i]+sigma_sumsoil[i-1]    
       
    return (sigma_sumsoil)
    

#Defining vectors -------------
class gamma_vectors:
    def __init__ (self):  
        self.s = []
        self.w = []
        
class Permeability_curve_Perm_vector:
    def __init__ (self):  
        self.K_tot = []
        self.K_min = []
        self.tau = []     

class Compression_vector:
    def __init__(self):
        self.Kappa = []
        self.Lambda = []
        self.Mu = []
        self.Tau_ref = []     
        self.OCR = []
                
class Water_retention_curve_vectors:
    def __init__ (self):  
        self.vG_WCR = []
        self.vG_WCS = []
        self.vG_a = []
        self.vG_alpha = []
        self.vG_m = []
        self.vG_mod = []
        self.vG_n = []

class Bc_top:
    def __init__ (self):    
        self.pressure = []
        
class Bc_bottom:
    def __init__ (self):   
        self.pressure = []  
        
class Soil_stresses:
    def __init__(self):
        self.sigma_total_node = []
        self.sigma_effective_node = []
        self.sigma_effective_node_ini = []
        
class Soil_states:
    def __init__(self):
        self.e_ini = []
        self.theta_ini= []
        self.e_out = []
        self.S_e = []
        self.Plastic_point =[]
        
                
class Soil_calc:
    def __init__(self): 
        self.de_dtheta_drying = []
        self.d2e_dtheta2_drying = []
        self.de_dtheta_wetting = []
        self.d2e_dtheta2_wetting = []
        self.theta = []

class swp:
    def __init__(self): 
        self.depth = []
        self.wcontent = []