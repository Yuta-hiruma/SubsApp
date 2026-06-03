Soil_parameters ={
    'aquifer':{      
        'Material_parameters':{
        #'gamma_s' : 2.65,
        'gamma_s' : 13.7,
        'gamma_w' : 1,
        'gamma_o' : 1.7,
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
        'kappa' : 0.00173718,   
        'lambda'  : 0.001736,
        'mu'  : 0.0001,
        'Tau_ref'  : 1,
        'Pg' : 200,
        'OCR' : 1.5,
        'c' : 10,
        'phi' : 20,
        },
        
        'Water_retention_curve':{
        'vG_WCS' : 0.35,
        'vG_WCR' : 0.01,
        'vG_alpha' : 0.030,
        'vG_n' : 2.88,
        'vG_m' : 0.18,
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
        'K_tot' : 0.0001,
        'K_min' : 0.0001,
        'tau' : 0.5
        }
        
        }
    }

Soil_parameters['confining_layer'] = {
        'Material_parameters':{
        #'gamma_s' : 2.65,
        'gamma_s' : 13.7,
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
        'kappa' : 0.0868,   
        'lambda'  : 0.1302,
        'mu'  : 0.0028229,
        'Tau_ref'  : 1,
        'Pg' : 200,
        'OCR' : 1.5,
        'c' : 0,
        'phi' : 30,
        },
        
        'Water_retention_curve':{
        'vG_WCS' : 0.5,
        'vG_WCR' : 0.00,
        'vG_alpha' : 0.0097,
        'vG_n' : 1.36,
        'vG_m' : (1-1/1.31),
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
        'K_tot' : 0.001,
        'K_min' : 0.0001,
        'tau' : 0.5
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