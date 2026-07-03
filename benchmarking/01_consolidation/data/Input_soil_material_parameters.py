Soil_parameters ={
    'aquifer':{      
        'Material_parameters':{
        #'gamma_s' : 2.65,
        #'gamma_s' : 13.7,
        'gamma_s' : 2.31,
        'gamma_w' : 1,
        },
        

        'Compression':{
        'kappa' : 0.00173718,   
        'lambda'  : 0.001736,
        'mu'  : 0.0001,
        'Tau_ref'  : 1,
        'Pg' : 400,
        'OCR' : 1.5,
        'c' : 10,
        'phi' : 20,
        },
        
        'Water_retention_curve':{
        'vG_WCS' : 0.53,
        'vG_WCR' : 0.00,
        'vG_alpha' : 0.023,
        'vG_n' : 1.28,
        'vG_m' : (1-1/1.28),
        'vG_a' : 2000,
        'vG_mod' : 1,
        'Ksatfit' : 87.4,
        'Lexp' : -4.493,
        'H_enpr' :-10.0,
        'Ksattexm' : 1.00,
        'Bdens': 1300.0
        },
        
        'Permeability_curve_Perm':{
        'K_tot' : 0.001,
        'K_min' : 0.0001,
        'tau' : 0.5
        }
        
        }
    }

Soil_parameters['confining_layer'] = {
        'Material_parameters':{
        #'gamma_s' : 2.65,
        #'gamma_s' : 13.7,
        'gamma_s' : 2.31,
        'gamma_w' : 1,
        },
        
        'Compression':{
        'kappa' : 0.0868,   
        'lambda'  : 0.1302,
        'mu'  : 0.0028229,
        'Tau_ref'  : 1,
        'Pg' : 400,
        'OCR' : 1.5,
        'c' : 0,
        'phi' : 30,
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
        'K_tot' : 0.2,
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