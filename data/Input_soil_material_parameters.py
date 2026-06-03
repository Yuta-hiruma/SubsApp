
Soil_parameters = {
    'Peat':{
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


        'Compression': {
        'kappa': 0.04369966075434302,
        'lambda': 0.09092858881212801,
        'mu': 0.0051214402063013975,
        'Tau_ref': 1,
        'Pg': 400,
        'OCR': 1.9350401988239156,
        'ccott': 70.58445188206727,
        },
        
        'Water_retention_curve':{
        'vG_WCS' : 0.89,
        'vG_WCR' : 0.00,
        'vG_alpha' : 0.0097,
        'vG_n' : 1.36,
        'vG_m' : (1-1/1.36),
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
        'K_tot' : 0.03373876734545294,
        'K_min' : 0.0001,
        'tau' : 0.5
        },

        }
    }

Soil_parameters['Peat1_SWP'] = {
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


        'Compression': {
        'kappa': 0.04369966075434302,
        'lambda': 0.09092858881212801,
        'mu': 0.0051214402063013975,
        'Tau_ref': 1,
        'Pg': 400,
        'OCR': 1.9350401988239156,
        'ccott': 70.58445188206727,
        },
        
        'Water_retention_curve':{
        'vG_WCS' : 0.849,
        'vG_WCR' : 0.01,
        'vG_alpha' : 0.0097,
        'vG_n' : 1.36,
        'vG_m' : (1-1/1.31),
        'vG_a' : 50000,
        'vG_mod' : 1,
        'Ksatfit' : 1.46,
        'Lexp' : -10.502,
        'H_enpr' :-10.0,
        'Ksattexm' : 3.4,
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
        'K_tot' : 0.03373876734545294,
        'K_min' : 0.0001,
        'tau' : 0.5
        },
 
        'Emission':{
        'T_min' : -10,          
        'a' : 0.05,  
        'alpha' : 2.59,           
        'beta' : 1.84,
        }        
        
        }

Soil_parameters['Peat2_SWP'] = {
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


        'Compression': {
        'kappa': 0.04369966075434302,
        'lambda': 0.09092858881212801,
        'mu': 0.0051214402063013975,
        'Tau_ref': 1,
        'Pg': 400,
        'OCR': 1.9350401988239156,
        'ccott': 70.58445188206727,
        },
        
        'Water_retention_curve':{
        'vG_WCS' : 0.58,
        'vG_WCR' : 0.01,
        'vG_alpha' : 0.0126,
        'vG_n' : 1.31,
        'vG_m' : (1-1/1.31),
        'vG_a' : 50000,
        'vG_mod' : 1,
        'Ksatfit' : 35.9,
        'Lexp' : -4.493,
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
        'K_tot' : 0.03373876734545294,
        'K_min' : 0.0001,
        'tau' : 0.5
        },
       
        
        }


 
Soil_parameters ['Clay_top_SWP']= {      
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

        'Compression': {
        'kappa': 0.04369966075434302,
        'lambda': 0.09092858881212801,
        'mu': 0.0051214402063013975,
        'Tau_ref': 1,
        'Pg': 400,
        'OCR': 1.9350401988239156,
        'ccott': 70.58445188206727,
        },
        
        'Water_retention_curve':{
        'vG_WCS' : 0.53,
        'vG_WCR' : 0.01,
        'vG_alpha' : 0.023,
        'vG_n' : 1.28,
        'vG_m' : (1-1/1.28),
        'vG_a' : 2000,
        'vG_mod' : 1,
        'Ksatfit' : 87.4,
        'Lexp' : -4.493,
        'H_enpr' :-10.0,
        'Ksattexm' : 3.83,
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
        'K_tot' : 0.03373876734545294,
        'K_min' : 0.0001,
        'tau' : 0.5
        },

        
        }
Soil_parameters ['Clay2_SWP']= {      
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

        'Compression': {
        'kappa': 0.04369966075434302,
        'lambda': 0.09092858881212801,
        'mu': 0.0051214402063013975,
        'Tau_ref': 1,
        'Pg': 400,
        'OCR': 1.9350401988239156,
        'ccott': 70.58445188206727,
        },
        
        'Water_retention_curve':{
        'vG_WCS' : 0.44,
        'vG_WCR' : 0.00,
        'vG_alpha' : 0.014,
        'vG_n' : 1.12,
        'vG_m' : 0.18,
        'vG_a' : 2000,
        'vG_mod' : 1,
        'Ksatfit' : 2.12,
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
        'K_tot' : 0.03373876734545294,
        'K_min' : 0.0001,
        'tau' : 0.5
        },

        
        }
Soil_parameters ['Sand_SWP']= {      
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
        'kappa' : 0.0868,   
        'lambda'  : 0.1302,
        'mu'  : 0.0028229,
        'Tau_ref'  : 1,
        'Pg' : 400,
        'OCR' : 10,
        'ccott' : 50
        },
        
        'Water_retention_curve':{
        'vG_WCS' : 0.33,
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
        'K_tot' : 0.03373876734545294,
        'K_min' : 0.0001,
        'tau' : 0.5
        },

        
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