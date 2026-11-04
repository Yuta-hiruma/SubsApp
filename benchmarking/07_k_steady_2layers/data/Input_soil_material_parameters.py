Soil_parameters ={
    'Peat':{
        'Material_parameters':{
         'gamma_s' : 2.31,
        'gamma_w' : 1,
        },
        
        'Compression': {
        'kappa': 0.2,
        'lambda': 0.40,
        'mu': 0.025,
        'Tau_ref': 1,
        'Pg': 400,
        'OCR': 1.3,
        'c' : 10,
        'phi' : 20
        },
        
        'Water_retention_curve':{
        'vG_WCS' : 0.89,
        'vG_WCR' : 0.00,
        'vG_alpha' : 0.0097,
        'vG_n' : 1.28,
        'vG_m' : (1-1/1.28),
        'vG_a' : 0,
        'vG_mod' : 1,
        'Ksatfit' : 1.46,
        'Lexp' : -10.502,
        'H_enpr' :-10.0,
        'Ksattexm' : 3.08,
        'Bdens': 1300.0
        },
        
        'Permeability_curve_Perm':{
        'K_tot' : 0.1,
        'K_min' : 0.0001,
        'tau' : 0.5
        }       
        
        }
        }

Soil_parameters ['Clay1']= {      
        'Material_parameters':{
        'gamma_s' : 2.31,
        'gamma_w' : 1,
        },

        'Compression': {
        'kappa' : 0.0868,   
        'lambda'  : 0.1302,
        'mu'  : 0.0028229,
        'Tau_ref': 1,
        'Pg': 400,
        'OCR': 1.3,
        'c' : 10,
        'phi' : 20
        },
        
        'Water_retention_curve':{
        'vG_WCS' : 0.53,
        'vG_WCR' : 0.00,
        'vG_alpha' : 0.023,
        'vG_n' : 1.28,
        'vG_m' : (1-1/1.28),
        'vG_a' : 0,
        'vG_mod' : 1,
        'Ksatfit' : 87.4,
        'Lexp' : -4.493,
        'H_enpr' :-10.0,
        'Ksattexm' : 1.00,
        'Bdens': 1300.0
        },
        
        'Permeability_curve_Perm':{
        'K_tot' : 0.01,
        'K_min' : 0.0001,
        'tau' : 0.5
        }
        
        }

Soil_parameters ['Clay2']= {      
        'Material_parameters':{
        'gamma_s' : 2.31,
        'gamma_w' : 1,
        },

        'Compression': {
        'kappa' : 0.0868,   
        'lambda'  : 0.1302,
        'mu'  : 0.0028229,
        'Tau_ref': 1,
        'Pg': 400,
        'OCR': 1.3,
        'c' : 10,
        'phi' : 20
        },
        
        'Water_retention_curve':{
        'vG_WCS' : 0.53,
        'vG_WCR' : 0.00,
        'vG_alpha' : 0.023,
        'vG_n' : 1.28,
        'vG_m' : (1-1/1.28),
        'vG_a' : 0,
        'vG_mod' : 1,
        'Ksatfit' : 87.4,
        'Lexp' : -4.493,
        'H_enpr' :-10.0,
        'Ksattexm' : 1.00,
        'Bdens': 1300.0
        },
        
        'Permeability_curve_Perm':{
        'K_tot' : 0.015,
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




