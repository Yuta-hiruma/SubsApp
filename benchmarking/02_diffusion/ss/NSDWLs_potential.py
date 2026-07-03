# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 11:36:27 2022

@author: tomdegast
"""

import numpy as np
from scipy.integrate import solve_ivp

from NSDWLs_van_Genuchten import *

"Gravitational potential"

def pot_g (Soil_dimension):
    #pot_g  = g(z-z0)
    #g      = 9.81 m/s2 
    #z_0    = Reference level [m] = NAP
    #z      = Heigth [m]
    
    z = Soil_dimension.node_vector
    
    g = 9.81 
    
    pot_g = g * z
    
    return (pot_g)    

"Matrix potential [suction, adsorption]"

def pot_psi (Soil_dimension, Soil_states, Water_retention_curve_vector):
    WRC = Water_retention_curve_vector
    #psi_in = Soil_states.psi_in
    theta = Soil_states.theta
    
    pot_psi = van_genuchten(WRC, theta)    
    return (pot_psi)
    
"Osmotic potential"

def pot_os (Soil_dimension,Soil_states):
    #psi_o  = -c * nu * alpha * R * T_k
    #c      = solute concentration [mmol/kg]
    #nu     = number of particles in solutuion over mlecules of solute ( nu = 1 non- ionizing) ( nu = number of ions per moecule) 
    #aplha  = osmotic coefficient
    #R      = Gas constant (8,31 J/mol)
    #T_k    = temperature in Kelvinn
    dm = Soil_dimension.internode_distance 
    e = Soil_states.e_out
    
    dz = dm*(1+e)
    return (pot_os)
    
"Hydrostatic potential"

def pot_hydrostat (Soil_dimension,Soil_states):
    dm = Soil_dimension.internode_distance 
    e = Soil_states.e_out
    
    dz = dm*(1+e)
    return (pot_hydrostat)
    
"Overburden potential"
    
def Flow_Potential_prime(Flow_Potential,Bctop,Bcbottom,Soil_dimension):
    psi_node = Flow_Potential.total_node
    
    Flow_Potential_prime = psi_node*0
    return Flow_Potential_prime


class Flow_Potential:
    def __init__(self):
        self.total = []
        self.g = []
        self.os = []
        self.hydrostat = []
        self.matrix = []
