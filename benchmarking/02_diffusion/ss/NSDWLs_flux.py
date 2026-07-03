import numpy as np
from scipy.integrate import solve_ivp

#  Import plotting functionality
import matplotlib.pyplot as plt
from NSDWLs_calc import *
# Soil_dimension

# Soil_parameters
# Soil_states
# Soil_stresses


def flux(t,Soil_states,Soil_stresses,Permeability_cruve_Perm_vector,Shrinkage_curve_vector, Soil_dimension, Bctop, Bcbottom):
    #  Calculate spacing between points
    T           =  Soil_states.theta_ini
    N           = len(Soil_states.theta_ini)
    
    def equations(t, T):
        depth_in    = Soil_dimension.depth_out
        node_size   = Soil_dimension.node_size
        e_in        = Soil_states.e_ini
        e_out       = Soil_states.e_out 
        theta       = T
        
        #  Boundary conditions set explicitly.  This is probably redundant
        #  as these numbers are set before the function is called and the
        #  code takes the time derivative at these points to be zero, but
        #  we'll ensure the proper boundary conditions anyway.
        T[0] = Bctop.theta
        T[-1] = Bcbottom.theta
        #  Initialize time derivative vector
        Tprime = np.zeros( (N) )
        #  We are being lazy and forcing the temperature of the boundaries
        #  not to change by setting the derivative to zero at those points
        Tprime[0] = 0.0
        Tprime[-1] = 0.0
        #  Implement the diffusion equation in the interior points using the
        #  central difference formula.  This can be vectorized for better
        #  performance, but I am just using a loop here.
        
        e_in = e_out
        
        e_out = update_e (Shrinkage_curve_vector,theta)
                
        depth_out = update_height (depth_in, e_in, e_out)
        node_size = Update_node_size(depth_out,e_out)
        
        K = update_K(Soil_states,Permeability_cruve_Perm_vector, Shrinkage_curve_vector)
       
        #print (depth_out[0])
        
        for i in range(1, N - 1):
            Tprime[i] = K[i] * (T[i+1] - 2 * T[i] + T[i-1]) / node_size[i-1]*node_size[i]
        return Tprime
            
    #  Solve the equation by calling solve_ivp and return the solution
    sol = solve_ivp(equations, [t[0], t[-1]], T, t_eval = t)
    return sol
    

def update_e (Scv,theta):
    A_sh = Scv.A_sh
    B_sh = Scv.B_sh
    C_sh = Scv.C_sh
    theta = theta
        
    e_out = Scv.A_sh*(1+(theta**Scv.C_sh)/(Scv.B_sh**Scv.C_sh))**(1/Scv.C_sh)
    return e_out

def update_height (depth_in, e_in, e_out):
    depth_out = depth_in*0
    depth_out[-1] = depth_in[-1]
    for i in range(len(depth_in)-1):
        Vs_temp = (1-e_in[-i-1]/(1+e_in[-i-1]))*(depth_in[-i-2]-depth_in[-i-1])
        depth_out[-i-2] = depth_out[-i-1] + Vs_temp/(1-e_out[-i-1]/(1+e_out[-i-1]))      
    return (depth_out)


def Update_node_size (depth_in,e_out):
    node_size_out = e_out*0
    for i in range(len(depth_in)-1):
        node_size_out [-i-1] = depth_in[-i-1] - depth_in[-i-2]  
    return(node_size_out)
  


