# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 11:36:27 2022

@author: tomdegast
"""
import numpy as np
import scipy as sci

from NSDWLs_calc import *
"eq-1"

#   dm = dz/(1+e)

#   m is the 1D material coordinate
#   z is the Cartesian (real) coordinate
#   e is the void ratio
 
#def dz_dm (Soil_dimension,Soil_states):
#    dz = Soil_dimension.internode_distance 
#    e = Soil_states.e
#    
#    dm = dz/(1+e)
#    return (dm)

#def dm_dz (Soil_dimension,Soil_states):
#    dm = Soil_dimension.internode_distance 
#    e = Soil_states.e_out
#    
#    dz = dm*(1+e)
#    return (dz)


"~~~Water transport (Richards Equation)~~~"

"eq-2"

#   Theta/dt = d^2/dz^2 [K(psi+z+omega)]

#   ϕ is water potential ! not used
#   z is the gravimetric component
#   omega is an overburdencomponent (the excess pore pressure in head)
#   psi is a suction component 
#   Theta is the volumetric water content Theta = V_w/V_t
#   V_w and V_t are the voume of water and total volume
#   t is time
#   K is the hydraulic conductivity

#def K_potential (Soil_calc):
#    K           = Soil_calc.K
#    #dpsi_dm     = Soil_calc.dpsi_dm
#    SF1         = Soil_calc.SF1
#    SF2         = Soil_calc.SF2
#    dtheta_dm   = Soil_calc.dtheta_dm
#    
#    ' TEMP '
#    K = np.append(K,K[-1])
#    ' TEMP '
#    
#    K_potential = K*(dtheta_dm)
#    
#    return (K_potential)

#def K_gradient (Soil_calc,Soil_dimension,Bcbottom,Bctop):
#    K_potential          = Soil_calc.K_potential
#    internode_distance   = Soil_dimension.internode_distance
#    #print(internode_distance)
#    K_gradient = np.zeros(len(internode_distance))
#    
#    for i in range(len(internode_distance)):
#        K_gradient [i] = (K_potential[i]-K_potential[i-1]/internode_distance[i])
#    
#    K_gradient[0]  = Bctop.devaporation
#    #K_gradient[-1] = Bcbottom.dflow
#  
#    return (K_gradient)



"eq-3"

#   Omega = de/dtheta intergral(z0-z) gamma dz
#   gamma = (theta+gamma_s)/(1+e)

#   e is the void ratio
#   theta is the volumetric water ratio,V_w / V_s
#   V_w / V_s are the voume of water and volume of solids
#   Omega is an overburdencomponent (the excess pore pressure in head)
#   z is the depth
#   z0 is the depth at the surface
#   gamma is the volumetric weight
#   gamma_s = volumetric weight solids

def calc_stress (Soil_dimension,gamma_vector,Soil_states,Water_retention_curve_vector,Soil_stresses):
    # input:  theta_n, state par, simsettings, spatial, dm_in, ZZ, k, m
    # Output : psi_n (van genuchten), e_n, dedtheta_n, d2edtheta2_n, dZ_in, k, m 
    #Initialise
    #   Soil_domain
    node_size_initial = get_internode_distance(Soil_dimension)
    
    #    gamma_vector
    gamma_s     = np.array(gamma_vector.s)
    gamma_o     = np.array(gamma_vector.o)
    gamma_o_ratio = np.array(gamma_vector.o_ratio)
    gamma_w     = np.array(gamma_vector.w)
    
    #   WRC
    WCR = np.array(Water_retention_curve_vector.vG_WCR)
    WCS = np.array(Water_retention_curve_vector.vG_WCS)
    
    # States
    theta = np.array(Soil_states.theta)
    e = np.array(Soil_states.e)
    
    # Stresses
    psi_in = Soil_stresses.psi_ini

    # step 1 calculate overburden
     
    # step 2 calculate total stress
    " calculate stress based on e and the ratio between organic and non-organic components on the nodes"
    S_sat = (theta -WCR) / (WCS -WCR)
    sigma_total_node = np.zeros(np.shape(gamma_s))
    sigma_total_node = 9.81*((1-gamma_o_ratio)*gamma_s*(1/(e+1))+(gamma_o_ratio)*gamma_o*(1/(e+1))+S_sat*gamma_w*(1-(1/(e+1))))
    
    "calculate sum of total stress from nodes and project on nodes and internodes"
    sigma_total_sum_internode = np.zeros(np.shape(gamma_s))
    sigma_total_sum_node = np.zeros(np.shape(gamma_s))
    sigma_total_sum_internode = np.append(sigma_total_sum_internode,0)

    for i in range(len(sigma_total_node)):
       sigma_total_sum_internode[i+1] = node_size_initial[i]*sigma_total_node[i]+sigma_total_sum_internode[i]    

    for i in range(len(sigma_total_sum_node)):
        sigma_total_sum_node[i] = 0.5 * sigma_total_sum_internode[i] + 0.5 * sigma_total_sum_internode[i+1]

    # step 3 calculate water pressure(s)
    sigma_pressure_sum_internode = np.zeros(np.shape(sigma_total_sum_internode))
    sigma_pressure_sum_node = np.zeros(np.shape(sigma_total_sum_node))
    
    for i in range(len(sigma_total_sum_node)):
        if S_sat[i] < 1:
            sigma_pressure_sum_internode[i+1] = 0
        if S_sat[i] == 1:
            sigma_pressure_sum_internode[i+1] = node_size_initial[i]*9.81 + sigma_pressure_sum_internode[i]
        sigma_pressure_sum_node[i] = 0.5*sigma_pressure_sum_internode[i]+0.5*sigma_pressure_sum_internode[i+1]   

    # step 4 calculate suction
    Soil_stresses.psi = van_genuchten_h (Water_retention_curve_vector, Soil_stresses, Soil_states)


    # step 5 calculate effective stress
    psi = np.array(Soil_stresses.psi)
    
    sigma_effective_node = (sigma_total_sum_node)-sigma_pressure_sum_node
    sigma_effective_internode = sigma_total_sum_internode - sigma_pressure_sum_internode
    
    for i in range(len(sigma_total_sum_node)-1):
        sigma_effective_internode[i+1] = sigma_effective_internode[i+1]+ 0.5*psi[i]+0.5*psi[i+1]

    Soil_stresses.sigma_total_node              = np.array(sigma_total_node)
    Soil_stresses.sigma_total_sum_node          = np.array(sigma_total_sum_node)
    Soil_stresses.sigma_total_sum_internode     = np.array(sigma_total_sum_internode)
    Soil_stresses.sigma_pressure_sum_internode  = np.array(sigma_pressure_sum_internode)
    Soil_stresses.sigma_pressure_sum_node       = np.array(sigma_pressure_sum_node)
    Soil_stresses.sigma_effective_node          = np.array(sigma_effective_node)
    Soil_stresses.sigma_effective_internode     = np.array(sigma_effective_internode)

    return (Soil_stresses)

def update_sigma_total (node_size,gamma_vector,Soil_states,Water_retention_curve_vector,Bctop_in):
    node_size = node_size
    
    gamma_s     = np.array(gamma_vector.s)
    gamma_o     = np.array(gamma_vector.o)
    gamma_o_ratio = np.array(gamma_vector.o_ratio)
    gamma_w     = np.array(gamma_vector.w)
    
    WCR = np.array(Water_retention_curve_vector.vG_WCR)
    WCS = np.array(Water_retention_curve_vector.vG_WCS)
    
    # States
    theta = np.array(Soil_states.theta)
    e = np.array(Soil_states.e)
    n = e/(1+e)
    water_top = Bctop_in.pressure
    if water_top < 0 : 0
    
    S_sat = np.clip((theta - WCR) / (WCS - WCR), 0, 1)
    sigma_total_node = np.zeros(np.shape(gamma_s))
    sigma_total_node = 9.81*((1-gamma_o_ratio)*gamma_s*(1/(e+1))+(gamma_o_ratio)*gamma_o*(1/(e+1))+S_sat*gamma_w*(1-(1/(e+1))))
    
    sigma_total = np.zeros(len(sigma_total_node))
    sigma_total_in = np.zeros(len(sigma_total_node)+1)
     
    for i in range(len(sigma_total_node)):
        sigma_total_in[i+1] = sigma_total_node[i]*node_size[i]+sigma_total_in[i]#+water_top
        sigma_total[i] = 0.5*sigma_total_in[i]+0.5*sigma_total_in[i+1]

    
    return (sigma_total)

"eq-4"

#   - paper dtheta/dt = d/dm {K*[dpsi/dtheta*dtheta/dm+(1+e)+(theta+gamma_s)dpsi/dtheta+d2e/dtheta2*intergral(z0-z)](theta+gamma_s)dm*dtheta/dm]}

#   dtheta/dt = d/dm {K*}[dpsi/dm + SF1 + SF2*dtheta/dm]
#   SF1 = (1+e) - (theta+gamma_s)*de/dtheta
#   SF2 = d2e/dtheta2 * int_0^m(theta+gamma_s)dm
    
#def SF1(Soil_states,Soil_calc,gamma_vector):
#    theta = Soil_calc.theta_internode
#    e = Soil_calc.e_internode
#    gamma_s = gamma_vector.s
#    de_dtheta = Soil_calc.de_dtheta
#    # equation
#    "check sign convention"
#    SF1 = (1+e) -(theta+gamma_s)*de_dtheta
#    return SF1

#def SF2(Soil_states,Soil_calc,gamma_vector,Soil_dimension):
#    d2e_dtheta2  = Soil_calc.d2e_dtheta2
#    theta       = Soil_states.theta
#    gamma_w         = gamma_vector.w
#    gamma_s         = gamma_vector.s
#    gamma_o         = gamma_vector.o
#    o_ratio         = gamma_vector.o_ratio
#    internode_distance = Soil_dimension.internode_distance
#    gamma_s_net = o_ratio*gamma_o+(1-o_ratio)*gamma_s
#    
#    local = np.zeros(len(internode_distance)+1)
#    intergration = np.zeros(len(internode_distance)+1)
#    
#    for i in range(len(internode_distance)):
#        local      [i]      = theta[i]*gamma_w[i]+gamma_s_net[i]
#        intergration[i+1] = intergration[i]+local[i]
#        
#    print(intergration[9])    
#    SF2 = d2e_dtheta2*intergration
#        
#    return SF2 

#   K = K_sat*K_rel*K_dess (eq-6,eq-7,eq-8)

#   theta is the volumetric water ratio,V_w / V_s
#   t is time
#   m is the 1D material coordinate
#   K* = K/(1+e)
#   K is the hydraulic conductivity
#   psi is a suction component 
#   theta is the volumetric water ratio
#   e is the void ratio
#   gamma_s = volumetric weight solids
#   z is the depth
#   z0 is the depth at the surface

#def de_dtheta_drying(Shrinkage_curve_vector, Soil_states):
#    #initialise parameters
#    A_sh = Shrinkage_curve_vector.A_sh
#    B_sh = Shrinkage_curve_vector.B_sh
#    C_sh = Shrinkage_curve_vector.C_sh
#    theta= Soil_states.theta_internode
#    #equation
#    de_dtheta_drying = A_sh*(theta**C_sh/B_sh**C_sh + 1)**(1/C_sh)*theta**C_sh/(theta*(theta**C_sh + B_sh**C_sh)) 
#    return de_dtheta_drying

#def d2e_dtheta2_drying(Shrinkage_curve_vector, Soil_states):
#    #initialise parameters
#    A_sh = Shrinkage_curve_vector.A_sh
#    B_sh = Shrinkage_curve_vector.B_sh
#    C_sh = Shrinkage_curve_vector.C_sh
#    theta= Soil_states.theta_internode
#    #equation
#    d2e_dtheta2_drying = A_sh*(theta**C_sh/B_sh**C_sh + 1)**(1/C_sh)*theta**C_sh*B_sh**C_sh*(C_sh - 1)/(theta**2*(theta**C_sh + B_sh**C_sh)**2)
#    return d2e_dtheta2_drying

#def de_dtheta_wetting(Shrinkage_curve_vector, Soil_states):
#    # nu_h = void ratio at maximum swelling
#    # nu_z = the minimum void ratio of previouse shrinkage phase
#    # ksi_h = water content at maximum swelling
#    
#    nu_h = np.array(Shrinkage_curve_vector.n_eta)
#    nu_z = np.array(Shrinkage_curve_vector.n_eta)
#    ksi_h = np.array(Shrinkage_curve_vector.n_eta)
#    
#    theta = Soil_states.theta_internode
#    
#    de_dtheta_wetting = 2*(nu_h - nu_z)*(-theta + ksi_h)/ksi_h**2
#    return de_dtheta_wetting

#def d2e_dtheta2_wetting(Shrinkage_curve_vector, Soil_states):
#    # nu_h = void ratio at maximum swelling
#    # nu_z = the minimum void ratio of previouse shrinkage phase
#    # ksi_h = water content at maximum swelling
#    
#    nu_h = np.array(Shrinkage_curve_vector.n_eta)
#    nu_z = np.array(Shrinkage_curve_vector.n_eta)
#    ksi_h = np.array(Shrinkage_curve_vector.n_eta)
#    
#    theta = Soil_states.theta_internode
#    
#    d2e_dtheta2_wetting = 2*(nu_h - nu_z)*(-theta + ksi_h)/ksi_h**2
#    return d2e_dtheta2_wetting

#def dpsi_dm (Soil_stresses, Soil_dimension,Bctop,Bcbottom):
#    psi = Soil_stresses.psi
#    node_distance = Soil_dimension.node_distance
#    psi_top = Bctop.psi_top
#    
#    dpsi_dm_int = np.zeros(len(psi)+1)
#    for i in range(len(node_distance)):
#        dpsi_dm_int[i+1] = (psi[i]-psi[i+1])/node_distance[i]
#        
#    dpsi_dm_int[0]  = (psi[0]-psi_top)/node_distance[0]
#    dpsi_dm_int[-1] = 0
#    
#    
#    return dpsi_dm_int
#
#def dtheta_dm (Soil_states, Soil_dimension,Bctop,Bcbottom):
#    theta = Soil_states.theta
#    node_distance = Soil_dimension.node_distance
#    
#    dtheta_dm_int = np.zeros(len(theta)+1)
#    for i in range(len(node_distance)):
#        dtheta_dm_int[i+1] = (theta[i]-theta[i+1])/node_distance[i]
##    dtheta_dm_int[0]  = 0
#    dtheta_dm_int[-1] = 0
#    
#    
#    return dtheta_dm_int
    
"eq-5" "solved at eq-7"
#   K_sat = 10**(A_hcc*theta-B_hcc)

#   A_hcc and B_hcc are material parameters determined from experimental data
#   theta is the volumetric water ratio

"eq-6" "solved at eq-7"
#   K_rel = S_l**(delta)

#   Sl is the degree of saturation
#   delta is a material parameter and is related to the pore size distribution

"eq-7" "solved at eq-7"
#   K_dess = 1 if h-z > d_des
#   K_dess = (1-eps_dess)+eps_dess*(ksi_dess(1-S_l)) if h-z <d_dess

#   Sl is the degree of saturation
#   eps_dess is a material parameter
#   ksi_dess is a material parameter
#   z is the same as the gravitational potential (positive upwards), in this case with the zero datum at the base of the sludge
#   h is the total height of the sludge column
#   d_dess is the depth from the surface where desiccation and other surface effects act

def update_K (theta,e,S_e,Permeability_curve_Perm_vector,WRC):
    #initialise parameters    
    theta   = theta
    e       = e
    n = e/(1+e)
    S_l     = S_e
    
    for i in range(len(theta)):
        if theta[i] <= WRC.vG_WCR[i]:
            theta[i] = WRC.vG_WCR[i]
        elif theta[i] >= WRC.vG_WCS[i]:
            theta[i] = WRC.vG_WCS[i]
   
    A_HCC   = Permeability_curve_Perm_vector.A_HCC
    B_HCC   = Permeability_curve_Perm_vector.B_HCC
    delta   = Permeability_curve_Perm_vector.delta
    d_dess  = Permeability_curve_Perm_vector.desslength
    eps_dess= Permeability_curve_Perm_vector.epsilon
    ksi_dess= Permeability_curve_Perm_vector.eta
    K_out= Permeability_curve_Perm_vector.K_tot
    
    #equation 5
    K_sat = 10**(A_HCC*theta - B_HCC)
   #print ("K_sat",K_sat)
   #print("A_HCC",A_HCC)
   #print("theta",theta)
   #print("B_HCC",B_HCC)
    #equation 6
    S= theta/n
    K_rel = S**(delta)
    #equation 7
    K_dess = (1+eps_dess)+eps_dess*np.exp(ksi_dess*(1-S_l))     
    
    K_dess = np.ones(len(K_dess))
    
    #equation
    K_out = K_sat*K_rel*K_dess   
    
   
    
    return (K_out)

def update_K_h (theta,e,S_e,Permeability_curve_Perm_vector,WRC):
    #initialise parameters    

    n = e/(1+e)
    #tau = Permeability_curve_Perm_vector.tau
    tau = 0.5
    
    for i in range(len(theta)):
        if theta[i] <= WRC.vG_WCR[i]:
            theta[i] = WRC.vG_WCR[i]
        elif theta[i] >= WRC.vG_WCS[i]:
            theta[i] = WRC.vG_WCS[i]
   

    K_sat = Permeability_curve_Perm_vector.K_tot
    vG_n = WRC.vG_n
    vG_m = WRC.vG_m
    

    K_out = K_sat*(S_e**tau)*(1-(1-S_e**(1/vG_m)**vG_m))**2

    return (K_out)

def update_K_star(Soil_states):
    #initialise parameters  
    K = Soil_states.K
    e = Soil_states.e
    
    #equation 
    K_star = K/(1+e)
    return (K_star)

"~~~Soil deformation~~~"


"eq-8 - e_drying"
#   e = A_sh((theta**C_sh)/(B_sh**C_sh)+1)**(1/C_sh)

#   e is the void ratio
#   Ash is the minimum void ratio
#   Bsh is a parameter defining the slope
#   Csh is a parameter defining the transition between the linear portion and the minimum void ratio


#def update_e_drying (Shrinkage_curve_vector,Soil_states):
#    #initialise parameters
#    A_sh = Shrinkage_curve_vector.A_sh
#    B_sh = Shrinkage_curve_vector.B_sh
#    C_sh = Shrinkage_curve_vector.C_sh
#    theta = Soil_states.theta
#    #equation
#    e_out = A_sh*(1+(theta**C_sh)/(B_sh**C_sh))**(1/C_sh)
#    return e_out


"eq- 9"
#   A_sh = A_sh^0 *(1-(1/C_10)*log(sigma'/sigma'_0))

#   Ash^0 is defined as the minimum void ratio under zero overburden conditions
#   sigma´ is the current stress
#   sigma´_0 is the initial stress
#   C_10 is a material parameter

#def update_A_sh(Shrinkage_curve_vector,Soil_stresses):
#    #initialise parameters
#    A_sh0 = Shrinkage_curve_vector.A_sh0
#    C_10 = Shrinkage_curve_vector.C_10
#    sigma_eff = Soil_stresses.sigma_effective_node
#    sigma_eff0 = Soil_stresses.sigma_effective_node_ini
#    
#    #equation    
#    A_sh = A_sh0 *(1-(1/C_10)*np.log(sigma_eff/sigma_eff0))
#    
#    return A_sh

#"eq-9 ?"
#def update_B_sh (Shrinkage_curve_vector,Soil_states):
#    #initialise parameters
#    A_sh = np.array(Shrinkage_curve_vector.Ash)
#    S_l = np.array(Soil_states.S_e)
#    
#    #equation   
#    B_sh = S_l*A_sh
#    return B_sh    

"eq -10 - e_wetting"
#   e = nu_h - b(theta-ksi_h)**2
#   b= (nu_h-nu_z)/ksi_h**2

#   e is the void ratio
#   theta is the volumetric water ratio
#   nu_h is void ratio at the maximum swelling
#   ksi_h is the water content at the maximum swelling
#   nu_z is the minimum void ratio of the previous shrinkage phase

#def update_e_wetting (Shrinkage_curve_vector,Soil_states):
#    #initialise parameters
#    nu_h = Shrinkage_curve_vector.nu_h
#    nu_z = Shrinkage_curve_vector.nu_z
#    ksi_h = Shrinkage_curve_vector.ksi_h
#    
#    theta = Soil_states.theta 
#    
#    #equation
#    e_out = nu_h - ((nu_h-nu_z)/(ksi_h**2))*(theta-ksi_h)**2 
#    return e_out

#def check_wetting_drying_internode(Soil_calc):
#    de_dtheta_drying        =Soil_calc.de_dtheta_drying
#    d2e_dtheta2_drying      =Soil_calc.d2e_dtheta2_drying
#    de_dtheta_wetting       =Soil_calc.de_dtheta_wetting
#    d2e_dtheta2_wetting     =Soil_calc.d2e_dtheta2_wetting
#    
#    isdrying                =Soil_calc.isdrying
#    
#    Soil_calc.de_dtheta     = isdrying*de_dtheta_drying     +-1*(isdrying-1)*de_dtheta_wetting
#    Soil_calc.d2e_dtheta2   = isdrying*d2e_dtheta2_drying   +-1*(isdrying-1)*d2e_dtheta2_wetting
#      
#    
#    return (Soil_calc)

#def check_wetting_drying_node(Soil_calc):
#    e_drying                =Soil_calc.e_drying
#    e_wetting               =Soil_calc.e_wetting
#    isdrying                =Soil_calc.isdrying
#    e_out              = isdrying*e_drying             +-1*(isdrying-1)*e_wetting
#       
#    
#    return (e_out)
    

"~~~Material properties~~~"

"eq -11"
#   S_e = C(s)*(1/((1+(alpha_WRC*psi)**n_WRC)**m_WRC))
#   C(s) = 1 - ((ln(1+psi.aplha_WRC))/ln(2))
#   S_e = (theta - WRC)/(WCS-WCR)

#   Se is the effective degree of saturation
#   WCR is the residual (volumetric) water content
#   WCS is the water content at full saturation
#   alpha_WRC is a fitting parameters
#   n_WRC is a fitting parameters
#   m_WRC is a fitting parameters

def update_S_e (Water_retention_curve_vector,theta):
    vG_WCR = Water_retention_curve_vector.vG_WCR
    vG_WCS = Water_retention_curve_vector.vG_WCS
    theta = theta
    
    S_e = (theta - vG_WCR) / (vG_WCS - vG_WCR)
    return S_e


def van_genuchten (Water_retention_curve_vector, Soil_stresses, Soil_states):
    #Soil_states
    theta = np.array (Soil_states.theta)

    #shrinking
    Se = np.array(update_S_e(Water_retention_curve_vector, theta))
    
    # Water_retention_curve_vector
    vG_a        = Water_retention_curve_vector.vG_a
    vG_alpha    = Water_retention_curve_vector.vG_alpha
    vG_n        = Water_retention_curve_vector.vG_n
    vG_m        = Water_retention_curve_vector.vG_m

    #iterertavely slove psi
    psi_out = np.array(Soil_stresses.psi_ini)   
    
    #print (psi_out)
    for i in range(len(Se)):     
        psi_high = vG_a[i]
        psi_low = 0.0
        psi_avg = 0.5*psi_high+0.5*psi_low
        tolerance = np.abs(psi_high-psi_low)
        if Se[i] >= 1.0: 
            Se[i] = 1.0
            psi_out[i] = 0.0
            tolerance = 0.0
        if Se[i] <= 0.0 :
            Se[i] = 0.0
            psi_out[i] = vG_a[i]
            tolerance = 0.0
        while tolerance > 0.001:
            Se_high = (1-(np.log(1+psi_high/vG_a[i]))/np.log(2))*((1/(1+(vG_alpha[i]*psi_high)**vG_n[i])**vG_m[i]))
            Se_avg = (1-(np.log(1+psi_avg/vG_a[i]))/np.log(2))*((1/(1+(vG_alpha[i]*psi_avg)**vG_n[i])**vG_m[i]))
            Se_low = (1-(np.log(1+psi_low/vG_a[i]))/np.log(2))*((1/(1+(vG_alpha[i]*psi_low)**vG_n[i])**vG_m[i]))
            if Se_high == Se[i] : 
                psi_out[i] = psi_high 
                tolerance = 0.0
            if Se_avg == Se[i] : 
                psi_out[i] = psi_avg 
                tolerance = 0.0
            if Se_low == Se[i] : 
                psi_out[i] = psi_low 
                tolerance = 0.0
            psi_out[i] = psi_avg
            tolerance = np.abs(psi_high-psi_low)
            print (i, tolerance,"psi", psi_low, psi_avg, psi_high,"Se_avg", Se_avg,"Se[i]",Se[i])
            if Se_avg < Se[i] : 
                psi_high = psi_avg
                psi_avg = 0.5*psi_high+0.5*psi_low
            if Se_avg > Se[i] : 
                psi_low = psi_avg
                psi_avg = 0.5*psi_high+0.5*psi_low        
        #print(i,Se[i],psi_out[i])
    #print(psi_out)

    return -psi_out    
#-------------------------------------------
def van_genuchten_h (Water_retention_curve_vector, Soil_stresses, Soil_states):
    #Soil_states
    theta = np.array (Soil_states.theta)

    #shrinking
    Se = np.array(update_S_e(Water_retention_curve_vector, theta))
    # Water_retention_curve_vector
    vG_a        = Water_retention_curve_vector.vG_a
    vG_alpha    = Water_retention_curve_vector.vG_alpha
    vG_n        = Water_retention_curve_vector.vG_n
    vG_m        = Water_retention_curve_vector.vG_m
    #iterertavely slove psi
    psi_out = np.array(Soil_stresses.psi_ini)   

    h = (Se**(-1/vG_m)-1)/vG_alpha # pressure head(cm)
    psi_out = h * 0.09807                   #convert to kPa
    return -psi_out   
#-------------------------------------------
####~~~~ update dimensions ~~~ ###


def update_height (Soil_dimension, Soil_states):
    internode_vector = Soil_dimension.internode_vector
    e = Soil_states.e
    e_ini = Soil_states.e_ini
    
    internode_vector_out = internode_vector*0
    internode_vector_out[-1] = internode_vector[-1]
    for i in range(len(internode_vector)-1):
        Vs_temp = (1-e_ini[-i-1]/(1+e_ini[-i-1]))*(internode_vector[-i-2]-internode_vector[-i-1])
        internode_vector_out[-i-2] = internode_vector_out[-i-1] + Vs_temp/(1-e[-i-1]/(1+e[-i-1]))          
    return (internode_vector_out)
    
def get_node_coordinates(Soil_dimension):
    internode_vector = Soil_dimension.internode_vector
    node_vector = np.zeros((len(internode_vector)-1))
    for i in range(len(node_vector)):
        node_vector[i] = 0.5*(internode_vector[i]+internode_vector[i+1])
    return (node_vector)
    
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


#def update_theta_internode_from_theta_node (Soil_states,Soil_dimension):
#    theta_node = Soil_states.theta
#    node_size = Soil_dimension.internode_distance
#       
#    theta_internode = np.zeros(len(theta_node)+1)
#    
#    theta_internode[0] = theta_node[0]#-0.5*(theta_node[1]-theta_node[0])
#    theta_internode[-1] = theta_node[-1]#-0.5*(theta_node[-2]-theta_node[-1])
#    for i in range(len(theta_node)-1):
#        theta_internode[i+1] =  (theta_node[i]*node_size[i]+theta_node[i+1]*node_size[i+1])/(node_size[i]+node_size[i+1])
#    
#    for i in range(len(theta_internode)):
#        if theta_internode[i] < 1e-8:
#            theta_internode[i] = 1e-8
#    
#    return (theta_internode)

#def update_materialproperty_internode_from_node (materialproperty,Soil_dimension):
#    mp = materialproperty #note local input
#    node_size = Soil_dimension.internode_distance
#       
#    prop_internode = np.zeros(len(mp)+1)
#    
#    prop_internode[0] = mp[0]
#    prop_internode[-1] = mp[-1]
#    for i in range(len(mp)-1):
#        prop_internode[i+1] =  (mp[i]*node_size[i]+mp[i+1]*node_size[i+1])/(node_size[i]+node_size[i+1])
#    return (prop_internode)

#def update_Shrinkage_curve_vector_to_internode (Shrinkage_curve_vector,Soil_dimension):
#    
#    store = Shrinkage_curve_vectors()
#    
#    store.A_sh	 = update_materialproperty_internode_from_node (Shrinkage_curve_vector.A_sh,Soil_dimension)
#    store.A_sh0  = update_materialproperty_internode_from_node (Shrinkage_curve_vector.A_sh0,Soil_dimension)
#    store.B_sh   = update_materialproperty_internode_from_node (Shrinkage_curve_vector.B_sh,Soil_dimension)
#    store.C_10   = update_materialproperty_internode_from_node (Shrinkage_curve_vector.C_10,Soil_dimension)
#    store.C_sh   = update_materialproperty_internode_from_node (Shrinkage_curve_vector.C_sh,Soil_dimension)
#    store.C_10   = update_materialproperty_internode_from_node (Shrinkage_curve_vector.C_10,Soil_dimension)
#    store.n_eta  = update_materialproperty_internode_from_node (Shrinkage_curve_vector.n_eta,Soil_dimension)
#    store.v      = update_materialproperty_internode_from_node (Shrinkage_curve_vector.v,Soil_dimension)
#    store.nu_z   = update_materialproperty_internode_from_node (Shrinkage_curve_vector.nu_z,Soil_dimension)
#    store.nu_h  = update_materialproperty_internode_from_node (Shrinkage_curve_vector.nu_h,Soil_dimension)
#    store.ksi_h      = update_materialproperty_internode_from_node (Shrinkage_curve_vector.ksi_h,Soil_dimension)
#    return (store)

#def update_Permeability_curve_Perm_vector_to_internode (Permeability_curve_Perm_vector,Soil_dimension):
#    
#    store = Permeability_cruve_Perm_vectors()
#    
#    store.A_HCC		 = update_materialproperty_internode_from_node (Permeability_curve_Perm_vector.A_HCC,Soil_dimension)	
#    store.B_HCC      = update_materialproperty_internode_from_node (Permeability_curve_Perm_vector.B_HCC,Soil_dimension)   
#    store.delta      = update_materialproperty_internode_from_node (Permeability_curve_Perm_vector.delta,Soil_dimension)  
#    store.dess       = update_materialproperty_internode_from_node (Permeability_curve_Perm_vector.dess,Soil_dimension)  
#    store.desslength = update_materialproperty_internode_from_node (Permeability_curve_Perm_vector.desslength,Soil_dimension)
#    store.epsilon    = update_materialproperty_internode_from_node (Permeability_curve_Perm_vector.epsilon,Soil_dimension)
#    store.eta        = update_materialproperty_internode_from_node (Permeability_curve_Perm_vector.eta,Soil_dimension)
#    store.rel        = update_materialproperty_internode_from_node (Permeability_curve_Perm_vector.rel,Soil_dimension) 
#    store.K_tot      = update_materialproperty_internode_from_node (Permeability_curve_Perm_vector.K_tot,Soil_dimension) 
#    return (store)

#def update_Water_retention_curve_vector_to_internode (Water_retention_curve_vector,Soil_dimension):
#    
#    store = Water_retention_curve_vectors()
#    
#    store.vG_WCR	     = update_materialproperty_internode_from_node (Water_retention_curve_vector.vG_WCR,Soil_dimension)
#    store.vG_WCS     = update_materialproperty_internode_from_node (Water_retention_curve_vector.vG_WCS,Soil_dimension)
#    store.vG_a       = update_materialproperty_internode_from_node (Water_retention_curve_vector.vG_a,Soil_dimension)
#    store.vG_aplha   = update_materialproperty_internode_from_node (Water_retention_curve_vector.vG_alpha,Soil_dimension)
#    store.vG_m       = update_materialproperty_internode_from_node (Water_retention_curve_vector.vG_m,Soil_dimension)
#    store.vG_mod     = update_materialproperty_internode_from_node (Water_retention_curve_vector.vG_mod,Soil_dimension)
#    store.n          = update_materialproperty_internode_from_node (Water_retention_curve_vector.vG_n,Soil_dimension)
#    return (store)

#def update_gamma_vector_to_internode (gamma_vector,Soil_dimension):
#    
#    store = gamma_vectors()
#    
#    store.o	        = update_materialproperty_internode_from_node (gamma_vector.o,Soil_dimension)
#    store.o_ratio   = update_materialproperty_internode_from_node (gamma_vector.o_ratio,Soil_dimension)
#    store.s         = update_materialproperty_internode_from_node (gamma_vector.s,Soil_dimension)
#    store.w         = update_materialproperty_internode_from_node (gamma_vector.w,Soil_dimension)
#
#    return (store)

#def update_e_vector_to_internode (Soil_calc,Soil_dimension):
#    
#    store = Soil_calc
#    
#    store.e_drying_internode	= update_materialproperty_internode_from_node (Soil_calc.e_drying,Soil_dimension)
#    store.e_wetting_internode   = update_materialproperty_internode_from_node (Soil_calc.e_wetting,Soil_dimension)
#    
#    isdrying = store.isdrying
#    
#    store.e_internode              = isdrying*store.e_drying_internode             +-1*(isdrying-1)*store.e_wetting_internode
#    
#    return (store)

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
    # Fourth set is Plastic point #
    sigma_eff_max = in_input[3*nr_elements:4*nr_elements]
    # Fith set is Temperature #
    T = in_input[4*nr_elements:5*nr_elements]    
    # sixth set is Emmision #
    emm = in_input[5*nr_elements:6*nr_elements] 
    
    return(Pressure, theta, e,sigma_eff_max,T,emm)

def df_dz_non_uniform(Parameter,Bctop,Bcbottom,Soil_dimension):
   
    Parameter_loc = Parameter
    
    Bctop_loc = Bctop.pressure
    Bcbottom_loc = Bcbottom.pressure
    
    node_size_loc = 0*Parameter_loc
    df_dz = 0*Parameter_loc

    #for i in range(len(Parameter_loc)):
    #    node_size_loc[i]= Soil_dimension.internode_vector[i]-Soil_dimension.internode_vector[i+1]
    
    node_size_loc = Soil_dimension.node_distance
        
    #start   
    for i in range(1,len(df_dz)-1):
        df_dz[i] = (Parameter_loc[i+1]-Parameter_loc[i-1])/(0.5*(node_size_loc[i-1]+2*node_size_loc[i]+node_size_loc[i+1]))
    
    
    df_dz[0]=(Parameter_loc[1]-Bctop_loc)/(0.5*(3*node_size_loc[0]+node_size_loc[1]))
    
    df_dz[-1]=(Bcbottom_loc-Parameter_loc[-2])/(0.5*(3*node_size_loc[-1]+node_size_loc[-2]))
    
    
    return (df_dz)


def d2f_dz2_non_uniform(Parameter,Bctop,Bcbottom,Soil_dimension):
    
    Parameter_loc = Parameter
    
    Bctop_loc = Bctop
    Bcbottom_loc = Bcbottom
    
    node_size_loc = 0*Parameter_loc
    d2f_dz2 = 0*Parameter_loc

    #for i in range(len(Parameter_loc)):
    #    node_size_loc[i]= Soil_dimension.internode_vector[i]-Soil_dimension.internode_vector[i+1]
    
    node_size_loc = Soil_dimension.node_distance
    
    #start
    for i in range(1,len(d2f_dz2)-1):
        d2f_dz2[i] = (((Parameter_loc[i+1]-Parameter_loc[i])/node_size_loc[i])-((Parameter_loc[i]-Parameter_loc[i-1])/node_size_loc[i-1]))/(0.5*node_size_loc[i]+0.5*node_size_loc[i-1])
    
    d2f_dz2[0] = (((Parameter_loc[1]-Parameter_loc[0])/node_size_loc[0])-((Parameter_loc[0]-Bctop_loc)/node_size_loc[0]))/(0.5*node_size_loc[0]+0.5*node_size_loc[0])

    d2f_dz2[-1]=(((Bcbottom_loc-Parameter_loc[-1])/node_size_loc[-1])-((Parameter_loc[-1]-Parameter_loc[-2])/node_size_loc[-2]))/(0.5*node_size_loc[-1]+0.5*node_size_loc[-2])
 
    return (d2f_dz2)

def d2f_dz2_non_uniform_h(Parameter,Bctop,Bcbottom,Soil_dimension):
    
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

def d2f_dz2_non_uniform_2(Parameter,Bctop,Bcbottom,node_size):
    
    Parameter_loc = Parameter
    
    Bctop_loc = Bctop
    Bcbottom_loc = Bcbottom
    
    node_size_loc = 0*Parameter_loc
    d2f_dz2 = 0*Parameter_loc

    #for i in range(len(Parameter_loc)):
    #    node_size_loc[i]= Soil_dimension.internode_vector[i]-Soil_dimension.internode_vector[i+1]
    
    node_size_loc = node_size
    
    #start
    for i in range(1,len(d2f_dz2)-1):
        d2f_dz2[i] = (((Parameter_loc[i+1]-Parameter_loc[i])/node_size_loc[i])-((Parameter_loc[i]-Parameter_loc[i-1])/node_size_loc[i-1]))/(0.5*node_size_loc[i]+0.5*node_size_loc[i-1])
    
    d2f_dz2[0] = (((Parameter_loc[1]-Parameter_loc[0])/node_size_loc[0])-((Parameter_loc[0]-Bctop_loc)/node_size_loc[0]))/(0.5*node_size_loc[0]+0.5*node_size_loc[0])

    d2f_dz2[-1]=(((Bcbottom_loc-Parameter_loc[-1])/node_size_loc[-1])-((Parameter_loc[-1]-Parameter_loc[-2])/node_size_loc[-2]))/(0.5*node_size_loc[-1]+0.5*node_size_loc[-2])
 
    return (d2f_dz2)


def d2f_dz2_non_uniform_loc(Parameter,Bctop,Bcbottom,Soil_dimension):
    
    Parameter_loc = Parameter
    
    Bctop_loc = Bctop
    Bcbottom_loc = Bcbottom
    
    node_size_loc = 0*Parameter_loc
    d2f_dz2 = 0*Parameter_loc

    #for i in range(len(Parameter_loc)):
    #    node_size_loc[i]= Soil_dimension.internode_vector[i]-Soil_dimension.internode_vector[i+1]
    
    node_size_loc = Soil_dimension
    
    #start
    
    d2f_dz2 = (((Bcbottom_loc-Parameter_loc)/node_size_loc)-((Parameter_loc-Bctop_loc)/node_size_loc))/(0.5*node_size_loc+0.5*node_size_loc)

    return (d2f_dz2)


#def derrivative_shrinkage (SCV, theta, psi_m):
#    theta = theta
#    de_dtheta = np.zeros(np.shape(theta))
#    A = SCV.A_sh
#    B = SCV.B_sh
#    C = SCV.C_sh
#    psi_m = -psi_m
#    de_dtheta =((B**(-C))*(theta**(C-1))*((B**(-C))*(theta**C)+1)**(1/C-1))
#    
#    for i in range(len(psi_m)):
#        if psi_m[i] <= 0:
#            de_dtheta[i] = 0.0
#    
#    
#    return de_dtheta

def deSSc_dt(Compression_vector, S_e, sigma_eff,dpsi_dt,Water_retention_curve_vector, PP):
    S_e = S_e
    dpsi_dt_in = dpsi_dt
    #vG_WCR = Water_retention_curve_vector.vG_WCR
    #vG_WCS = Water_retention_curve_vector.vG_WCS
    
    Labmda = Compression_vector.Lambda
    kappa = Compression_vector.Kappa
    mu = Compression_vector.Mu
    PP = abs(PP)
    
    Tau_ref = Compression_vector.Tau_ref
    
    sigma_eff = sigma_eff
    dsigma_eff_dt = -dpsi_dt_in*S_e
    
    ccott = Compression_vector.ccott
    
    #print (PP[0],sigma_eff[0])
    
    Tau = (PP/sigma_eff)**((Labmda-kappa)/mu)
    
    #print (Tau)
    
    dElastic = dsigma_eff_dt*kappa/(ccott+sigma_eff)
    dPlastic = (1/Tau)*(mu/Tau_ref)*(sigma_eff/PP)**((Labmda-kappa)/mu)
    #dPlastic = (1/1)*(mu/1)*(sigma_eff/PP)**((Labmda-kappa)/mu)
    
    deps_dt = -(dElastic + dPlastic)

    dPP_dt = PP*(dPlastic/(Labmda-kappa))    
    
    #print("theta","%.1f" %  theta[3])
    #print("PP","%.1f" %  PP[3])
    #print("dpsi_dt","%.1f" %  dpsi_dt[3])
    #print("sigma_eff","%.1f" %  sigma_eff[3])
    #print("dsigma_eff_dt","%.1f" %  dsigma_eff_dt[3])
    #print("Tau","%.1f" %  Tau[3])
    #print("dElastic","%.1f" % dElastic[3])
    #print("dPlastic","%.1f" % dPlastic[3])
    #print("dPP","%.1f" %  dPP_dt[3])
    #rint("de_dt","%.1f" %  de_dt[3])
    return deps_dt, dPP_dt

def deSSc_dt_h(Compression_vector, Soil_states, S_e, sigma_eff,dsigma_eff_dt, dpsi_dt,Water_retention_curve_vector, PP):
    S_e = S_e
    dpsi_d_in = dpsi_dt
    e = Soil_states.e
    
    Lambda = Compression_vector.Lambda
    kappa = Compression_vector.Kappa
    mu = Compression_vector.Mu
    Lambda_m = Lambda/(1+e)
    kappa_m = kappa/(1+e)
    mu_m = mu/(1+e)


    Tau_ref = Compression_vector.Tau_ref
    OCR = Compression_vector.OCR
    sigma_eff = sigma_eff
    dsigma_eff_dt = -dpsi_dt*S_e
    ccott = Compression_vector.ccott

    exponent = (Lambda_m - kappa_m) / mu_m 
    dElastic = dsigma_eff_dt*kappa_m/(ccott+sigma_eff)
    dPlastic = 0
    #dPlastic = mu_m * (1 / Tau_ref) * (OCR) ** (exponent)

    deps_dt = -(dElastic + dPlastic)
    dPP_dt = PP*(dPlastic/Lambda_m - kappa_m) 
    #dPP_dt = np.clip(dPP_dt, -1e8, 1e8)    
    dn_dt = deps_dt/(1+e)**2

    return deps_dt, dn_dt, dPP_dt



def deSSc_dt_loc(Compression_vector, sigma_eff,dsigma_eff_dt,Water_retention_curve_vector, PP):
    #vG_WCR = Water_retention_curve_vector.vG_WCR
    #vG_WCS = Water_retention_curve_vector.vG_WCS
    
    Labmda = Compression_vector.Lambda
    kappa = Compression_vector.Kappa
    mu = Compression_vector.Mu
    ur = Compression_vector.unloadreload
    PP = abs(PP)
    
    Tau_ref = Compression_vector.Tau_ref
    
    sigma_eff = abs(sigma_eff)
    dsigma_loc_eff_dt = -(dsigma_eff_dt)
    
    ccott = Compression_vector.ccott
    
    #print (PP[0],sigma_eff[0])
    
    Tau = (PP/sigma_eff)**((Labmda-kappa)/mu)
    
    #print (Tau)
    if dsigma_loc_eff_dt < 0 : kappa = kappa*ur
    
    dElastic = dsigma_loc_eff_dt*kappa/(ccott+sigma_eff)
    dPlastic = (1/Tau)*(mu/Tau_ref)*(sigma_eff/PP)**((Labmda-kappa)/mu)
    #dPlastic = (1/1)*(mu/1)*(sigma_eff/PP)**((Labmda-kappa)/mu)
    
    deps_dt = -(dElastic + dPlastic)

    dPP_dt = PP*(dPlastic/(Labmda-kappa))    
    
    #print("theta","%.1f" %  theta[3])
    #print("PP","%.1f" %  PP[3])
    #print("dpsi_dt","%.1f" %  dpsi_dt[3])
    #print("sigma_eff","%.1f" %  sigma_eff[3])
    #print("dsigma_eff_dt","%.1f" %  dsigma_eff_dt[3])
    #print("Tau","%.1f" %  Tau[3])
    #print("dElastic","%.1f" % dElastic[3])
    #print("dPlastic","%.1f" % dPlastic[3])
    #print("dPP","%.1f" %  dPP_dt[3])
    #rint("de_dt","%.1f" %  de_dt[3])
    return deps_dt, dPP_dt

def deSSc_dt2(Compression_vector, sigma_eff,dsigma_eff_dt,Water_retention_curve_vector, PP):
    #vG_WCR = Water_retention_curve_vector.vG_WCR
    #vG_WCS = Water_retention_curve_vector.vG_WCS
    
    Labmda = Compression_vector.Lambda
    kappa = Compression_vector.Kappa
    mu = Compression_vector.Mu
    ur = Compression_vector.unloadreload
    PP = abs(PP)
    
    Tau_ref = Compression_vector.Tau_ref
    
    sigma_eff = sigma_eff
    dsigma_loc_eff_dt = -(dsigma_eff_dt)
    
    ccott = Compression_vector.ccott
    
    #print (PP[0],sigma_eff[0])
    
    Tau = (PP/sigma_eff)**((Labmda-kappa)/mu)
    
    #print (Tau)
    for i in range(len(dsigma_loc_eff_dt)):
        if dsigma_loc_eff_dt[i] < 0 : kappa[i] = kappa[i]*ur[i]
    
    dElastic = dsigma_loc_eff_dt*kappa/(ccott+sigma_eff)
    dPlastic = (1/Tau)*(mu/Tau_ref)*(sigma_eff/PP)**((Labmda-kappa)/mu)
    #dPlastic = (1/1)*(mu/1)*(sigma_eff/PP)**((Labmda-kappa)/mu)
    
    deps_dt = -(dElastic + dPlastic)

    dPP_dt = PP*(dPlastic/(Labmda-kappa))    
    
    #print("theta","%.1f" %  theta[3])
    #print("PP","%.1f" %  PP[3])
    #print("dpsi_dt","%.1f" %  dpsi_dt[3])
    #print("sigma_eff","%.1f" %  sigma_eff[3])
    #print("dsigma_eff_dt","%.1f" %  dsigma_eff_dt[3])
    #print("Tau","%.1f" %  Tau[3])
    #print("dElastic","%.1f" % dElastic[3])
    #print("dPlastic","%.1f" % dPlastic[3])
    #print("dPP","%.1f" %  dPP_dt[3])
    #rint("de_dt","%.1f" %  de_dt[3])
    return deps_dt, dPP_dt

def respiration_rate_curve (S_e,e,theta,dSatdt , dTemperature, Temperature, Water_retention_curve_vector,Emission_vector):
    T = Temperature
    T_min = Emission_vector.T_min
    dTemp = dTemperature
   
    a = Emission_vector.a # (Ratkowsky et al., 1983; Lloyd and Taylor, 1994; Bååth, 2018)
    
    dSedt = dSatdt
    Sat = S_e
    alpha = Emission_vector.alpha
    beta = Emission_vector.beta
    B_ab = sci.special.gamma(alpha+beta)/(sci.special.gamma(beta)*sci.special.gamma(alpha))   
    
    for i in range(len(Sat)):
        if Sat[i]>= 1 :
            Sat[i] = 1


    
    
    A_WFPS = sci.stats.beta.pdf(Sat,alpha,beta)
    A_T = (a*(T-T_min))**2
    Au = A_WFPS*A_T
    
    dA_WFPS =   B_ab * ((((Sat**(alpha-1))*(alpha-1)*((1-Sat)**(beta-1)))/Sat)-(((Sat**(alpha-1))*((1-Sat)**(beta-1))*(beta-1))/(1-Sat)))
    
    
    for i in range(len(Sat)):
        if Sat[i] >= 1:
            dA_WFPS[i] =0
        if Sat[i] <= 0:
            dA_WFPS[i] =0
    
    dA_WFPSdt =  dSedt*dA_WFPS
    
    dA_T = 2*(a**2)*(T-T_min)
    
    dA_Tdt = 2*(a**2)*(T-T_min)*dTemp
    
    for i in range(len(T)):
        if T[i] <= T_min[i]:
            dA_Tdt[i] = 0
            A_T[i] = 0
        
    dAudt = dA_WFPSdt*A_T + A_WFPS*dA_Tdt


    return Au, dAudt

def dtheta_de (theta,e):
    
    return dtheta_de

#~~~~~~~~~~~~~~classes used~~~~~~~~~~~~~~#

class gamma_vectors:
    def __init__ (self):  
        self.s = []
        self.w = []
        self.o = []
        self.o_ratio = []
        
class Permeability_cruve_Perm_vectors:
    def __init__ (self):  
        self.A_HCC = []
        self.B_HCC = []
        self.delta = []
        self.dess = []
        self.desslength = []
        self.epsilon = []
        self.eta = []
        self.rel = []
        
                
class Shrinkage_curve_vectors:
    def __init__ (self):  
        self.A_sh = []
        self.A_sh0 = []
        self.B_sh = []
        self.C_sh = []
        self.n_eta = []
        self.v = []
        self.C_10 = []
         
                
class Water_retention_curve_vectors:
    def __init__ (self):  
        self.vG_WCR = []
        self.vG_WCS = []
        self.vG_a = []
        self.vG_alpha = []
        self.vG_m = []
        self.vG_mod = []
        self.vG_n = []

class Emission_vector:
    def __init__ (self):  
        self.T_min = []
        self.a = []
        self.alpha = []
        self.beta = []

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
        self.time = []
        self.A_HCC = []
        self.B_HCC = []
        self.delta = []
        self.dess = []
        self.desslength = []
        self.epsilon = []
        self.eta = []
        self.rel = []
        
                
class Results_Shrinkage_curve_vectors:
    def __init__ (self):  
        self.time = []
        self.A_sh = []
        self.A_sh0 = []
        self.B_sh = []
        self.C_sh = []
        self.n_eta = []
        self.v = []
        self.C_10 = []
         
                
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
