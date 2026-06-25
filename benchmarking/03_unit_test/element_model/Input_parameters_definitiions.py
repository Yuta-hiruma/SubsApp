"""
"material parameters"
gamma_s -> [g/cm3]
gamma_w -> [g/cm3]
"""

"""
"shrinkage curve"
A_sh    -> [-]
B_sh    -> [-]
C_sh    -> [-]
n_eta   -> [-]
v       -> [-]
"""

"""
"Soil water retention curve"
vG_WCS     -> [-]
vG_WCR     -> [-]
vG_aplha   -> [-]
vG_n       -> [-]
vG_m       -> [-]
vG_a       -> [-]
vG_mod     -> [-]
"""

"""
"permeability curves" 
Perm_a              -> [-]
Perm_b              -> [-]
Perm_rel            -> true (1) or false (0) tag for unsaturated permeability
Perm_delta          ->  parameter for unsaturated permeability
Perm_dess           ->  true (1) or false (0) tag for surface enhanced permeability
Perm_epsilon        -> parametWCer for surface enhanced permeability
Perm_eta            -> parameter for surface enhanced permeability
Perm_desslength     -> [cm] parameter for surface enhanced permeability - depth dessication occurs until
"""

"""
"settlement due to overburden"
overburden  -> true (1) or false (0) tag for overburden compression
C10         -> [-]
sigma_init  -> kPa
"""

"""
"material parameters"
gamma_s = 2.31
gamma_w = 1
"""

"""
"shrinkage curve"
A_sh    = 0.48
B_sh    = 0.48
C_sh    = 4.47
n_eta   = 2
v       = 2
"""

"""
"Soil water retention curve"
vG_WCS     = 5.91
vG_WCR     = 0.00
vG_aplha   = 0.92
vG_n       = 1.15
vG_m       = 1-(1/vG_n)
vG_a       = 10000
vG_mod     = 1
"""

"""
"permeability curves" 
Perm_a              = 0.783
Perm_b              = 5
Perm_rel            = 1
Perm_delta          = 3
Perm_dess           = 1
Perm_epsilon        = 0.05
Perm_eta            = 5
Perm_desslength     = 10
"""

Soil_parameters = {
    'Peat_sample':{
        'Material_parameters':{
         'gamma_s' : 2.31,
        'gamma_w' : 1,
        'gamma_o' : 1.17,
        'gamma_o_ratio' : 1.0
        },
        
        'Shrinkage_curve':{
        'A_sh0'  : 1,     
        'A_sh'  : 1,
        'B_sh'  : 1,
        'C_sh'  : 1,
        'n_eta' : 1,
        'v'     : 1,
        'C_10'  : 10,
        'nu_z'  :1,
        'nu_h'  :1,
        'ksi_h' :1,

        },


        'Compression':{
        'kappa': 0.087,
        'lambda': 0.13,
        'mu': 0.005,
        'Tau_ref'  : 1,
        'Pg' : 400,
        'OCR' : 1.5,
        'c' : 1,
        'phi' : 5,
        'K0' : 0.1,
        'nu' : 0.1,
        },
        
        'Water_retention_curve':{
        'vG_WCS' : 0.89,
        'vG_WCR' : 0.00,
        'vG_alpha' : 0.11,
        'vG_n' : 1.23,
        'vG_m' : (1-1/1.23),
        'vG_a' : 50000,
        'vG_mod' : 1,
        'Ksatfit' : 1.46,
        'Lexp' : -10.502,
        'H_enpr' :-10.0,
        'Ksattexm' : 3.08,
        'Bdens': 1300.0
        },
        
        'Permeability_curve_Perm':{
        'A_HCC' : 1.6,          
        'B_HCC' : 4.4,  
        'rel' : 1,           
        'delta' : 4,         
        'dess' : 1,         
        'epsilon' : 1,       
        'eta' : 1,           
        'desslength' : 1,
        'K_tot' : 0.1,
        'tau' : 0.5
        },
 
        'Emission':{
        'T_min' : -10,          
        'a' : 0.05,  
        'alpha' : 2.59,           
        'beta' : 1.84,
        }

        }
    }

Soil_parameters['silty_clay_loam'] = {
        'Material_parameters':{
         'gamma_s' : 2.59,
        'gamma_w' : 1,
        'gamma_o' : 1.17,
        'gamma_o_ratio' : 0.5
        },
        
        'Shrinkage_curve':{
        'A_sh0'  : 1,     
        'A_sh'  : 1,
        'B_sh'  : 1,
        'C_sh'  : 1,
        'n_eta' : 1,
        'v'     : 1,
        'C_10'  : 10,
        'nu_z'  :1,
        'nu_h'  :1,
        'ksi_h' :1,

        },


        'Compression':{
        'kappa': 0.087,
        'lambda': 0.13,
        'mu': 0.005,
        'Tau_ref'  : 1,
        'Pg' : 400,
        'OCR' : 1.5,
        'c' : 1,
        'phi' : 5,
        'K0' : 0.1,
        'nu' : 0.1,
        },
        
        'Water_retention_curve':{
        'vG_WCS' : 0.68,
        'vG_WCR' : 0.01,
        'vG_alpha' : 0.11,
        'vG_n' : 1.23,
        'vG_m' : (1-1/1.23),
        'vG_a' : 50000,
        'vG_mod' : 1,
        'Ksatfit' : 1.46,
        'Lexp' : -10.502,
        'H_enpr' :-10.0,
        'Ksattexm' : 3.08,
        'Bdens': 1300.0
        },
        
        'Permeability_curve_Perm':{
        'A_HCC' : 1.6,          
        'B_HCC' : 4.4,  
        'rel' : 1,           
        'delta' : 4,         
        'dess' : 1,         
        'epsilon' : 1,       
        'eta' : 1,           
        'desslength' : 1,
        'K_tot' : 0.1,
        'tau' : 0.5
        },
 
        'Emission':{
        'T_min' : -10,          
        'a' : 0.05,  
        'alpha' : 2.59,           
        'beta' : 1.84,
        }        
        
        }

Soil_parameters['silty clay'] = {
        'Material_parameters':{
         'gamma_s' : 2.577,
        'gamma_w' : 1,
        'gamma_o' : 1.17,
        'gamma_o_ratio' : 0.5
        },
        'Shrinkage_curve':{
        'A_sh0'  : 1,     
        'A_sh'  : 1,
        'B_sh'  : 1,
        'C_sh'  : 1,
        'n_eta' : 1,
        'v'     : 1,
        'C_10'  : 10,
        'nu_z'  :1,
        'nu_h'  :1,
        'ksi_h' :1,

        },


        'Compression':{
        'kappa': 0.087,
        'lambda': 0.13,
        'mu': 0.005,
        'Tau_ref'  : 1,
        'Pg' : 400,
        'OCR' : 1.5,
        'c' : 1,
        'phi' : 5,
        'K0' : 0.1,
        'nu' : 0.1,
        },
        
        'Water_retention_curve':{
        'vG_WCS' : 0.75,
        'vG_WCR' : 0.01,
        'vG_alpha' : 0.11,
        'vG_n' : 1.23,
        'vG_m' : (1-1/1.23),
        'vG_a' : 50000,
        'vG_mod' : 1,
        'Ksatfit' : 1.46,
        'Lexp' : -10.502,
        'H_enpr' :-10.0,
        'Ksattexm' : 3.08,
        'Bdens': 1300.0
        },
        
        'Permeability_curve_Perm':{
        'A_HCC' : 1.6,          
        'B_HCC' : 4.4,  
        'rel' : 1,           
        'delta' : 4,         
        'dess' : 1,         
        'epsilon' : 1,       
        'eta' : 1,           
        'desslength' : 1,
        'K_tot' : 0.1,
        'tau' : 0.5
        },
 
        'Emission':{
        'T_min' : -10,          
        'a' : 0.05,  
        'alpha' : 2.59,           
        'beta' : 1.84,
        }        
        
        }

Soil_parameters['unknown'] = {
        'Material_parameters':{
         'gamma_s' : 2.56,
        'gamma_w' : 1,
        'gamma_o' : 1.17,
        'gamma_o_ratio' : 0.5
        },
        'Shrinkage_curve':{
        'A_sh0'  : 1,     
        'A_sh'  : 1,
        'B_sh'  : 1,
        'C_sh'  : 1,
        'n_eta' : 1,
        'v'     : 1,
        'C_10'  : 10,
        'nu_z'  :1,
        'nu_h'  :1,
        'ksi_h' :1,

        },


        'Compression':{
        'kappa': 0.087,
        'lambda': 0.13,
        'mu': 0.005,
        'Tau_ref'  : 1,
        'Pg' : 400,
        'OCR' : 1.5,
        'c' : 1,
        'phi' : 5,
        'K0' : 0.1,
        'nu' : 0.1,
        },
        
        'Water_retention_curve':{
        'vG_WCS' : 0.75,
        'vG_WCR' : 0.01,
        'vG_alpha' : 0.11,
        'vG_n' : 1.23,
        'vG_m' : (1-1/1.23),
        'vG_a' : 50000,
        'vG_mod' : 1,
        'Ksatfit' : 1.46,
        'Lexp' : -10.502,
        'H_enpr' :-10.0,
        'Ksattexm' : 3.08,
        'Bdens': 1300.0
        },
        
        'Permeability_curve_Perm':{
        'A_HCC' : 1.6,          
        'B_HCC' : 4.4,  
        'rel' : 1,           
        'delta' : 4,         
        'dess' : 1,         
        'epsilon' : 1,       
        'eta' : 1,           
        'desslength' : 1,
        'K_tot' : 0.1,
        'tau' : 0.5
        },
 
        'Emission':{
        'T_min' : -10,          
        'a' : 0.05,  
        'alpha' : 2.59,           
        'beta' : 1.84,
        }        
        
        }

Soil_parameters ['sand']= {      
        'Material_parameters':{
        'gamma_s' : 2.31,
        'gamma_w' : 1,
        'gamma_o' : 1.5,
        'gamma_o_ratio' : 0.0
        },
        
        'Shrinkage_curve':{
        'A_sh0' : 0.48,   
        'A_sh'  : 0.9,
        'B_sh'  : 0.79,
        'C_sh'  : 2.85,
        'n_eta' : 1,
        'v'     : 1,
        'C_10'  : 10,
        'nu_z'  :1,
        'nu_h'  :1,
        'ksi_h' :1,
        },

        'Compression':{
        'kappa': 0.087,
        'lambda': 0.13,
        'mu': 0.005,
        'Tau_ref'  : 1,
        'Pg' : 400,
        'OCR' : 1.5,
        'c' : 1,
        'phi' : 5,
        'K0' : 0.1,
        'nu' : 0.1,
        },
        
        'Water_retention_curve':{
        'vG_WCS' : 0.33,
        'vG_WCR' : 0.01,
        'vG_alpha' : 0.11,
        'vG_n' : 1.23,
        'vG_m' : (1-1/1.23),
        'vG_a' : 2000,
        'vG_mod' : 1,
        'Ksatfit' : 17.4,
        'Lexp' : -4.493,
        'H_enpr' :-10.0,
        'Ksattexm' : 1.00,
        'Bdens': 1300.0
        },
        
        'Permeability_curve_Perm':{
        'A_HCC' : 0.783,          
        'B_HCC' : 5,  
        'rel' : 1,           
        'delta' : 3,         
        'dess' : 1,         
        'epsilon' : 0.05,       
        'eta' : 5,           
        'desslength' : 110,
        'K_tot' : 0.01,
        'tau' : 0.5
        },

        'Emission':{
        'T_min' : -10,          
        'a' : 0.00,  
        'alpha' : 2.59,           
        'beta' : 1.84,
        } 
        
        }



    
        
class Soil_domain:
    def __init__ (self):
        self.Depth_top = []
        self.Depth_bottom = []
        self.precision = []
        self.node_initial = []




"""
"settlement due to overburden"
overburden  = 1
C10         = 100
sigma_init  = 0.01
"""