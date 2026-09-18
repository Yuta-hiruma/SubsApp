Soil_parameters ={
    'Clay':{
        'Material_parameters':{
        'gamma_s' : 2.32,
        'gamma_w' : 1,
        },

        'Compression': {
        'kappa': 0.15,
        'lambda': 0.25,
        'mu': 0.007,
        'Tau_ref': 1,
        'OCR': 1.0,
        'c' : 10,
        'phi' : 20
        },
        
        'Water_retention_curve':{
        'vG_WCS' : 0.66645,
        'vG_WCR' : 0.00,
        'vG_alpha' : 0.05,
        'vG_n' : 1.1,
        'vG_m' : (1-1/1.1),
        'vG_a' : 10000,
        'vG_mod' : 1,
        'Ksatfit' : 1.46,
        'Lexp' : -10.502,
        'H_enpr' :-10.0,
        'Ksattexm' : 3.08,
        'Bdens': 1300.0
        },
        
        'Permeability_curve_Perm':{
        'K_tot' : 0.2,
        'K_min' : 0.000001,
        'tau' : 1
        }       
        
        }
        }
    
        
class Soil_domain:
    def __init__ (self):
        self.Depth_top = []
        self.Depth_bottom = []
        self.precision = []
        self.node_initial = []




