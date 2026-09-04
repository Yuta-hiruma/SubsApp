Soil_parameters ={
    'Clay':{
        'Material_parameters':{
        'gamma_s' : 2.32,
        'gamma_w' : 1,
        },

        'Compression': {
        'kappa': 0.26,
        'lambda': 0.42,
        'mu': 0.05,
        'Tau_ref': 1,
        'OCR': 1.0,
        'c' : 2,
        'phi' : 10
        },
        
        'Water_retention_curve':{
        'vG_WCS' : 0.735,
        'vG_WCR' : 0.00,
        'vG_alpha' : 0.05,
        'vG_n' : 1.15,
        'vG_m' : (1-1/1.15),
        'vG_a' : 50000,
        'vG_mod' : 1,
        'Ksatfit' : 1.46,
        'Lexp' : -10.502,
        'H_enpr' :-10.0,
        'Ksattexm' : 3.08,
        'Bdens': 1300.0
        },
        
        'Permeability_curve_Perm':{
        'K_tot' : 0.001,
        'K_min' : 0.00001,
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




