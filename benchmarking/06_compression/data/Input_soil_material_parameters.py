Soil_parameters ={
    'Clay':{
        'Material_parameters':{
        'gamma_s' : 2.31,
        'gamma_w' : 1,
        },

        'Compression': {
        'kappa': 0.01,
        'lambda': 0.28199239,
        'mu': 0.07,
        'Tau_ref': 1,
        'OCR': 1.3,
        'c' : 10,
        'phi' : 20
        },
        
        'Water_retention_curve':{
        'vG_WCS' : 0.735,
        'vG_WCR' : 0.00,
        'vG_alpha' : 0.023,
        'vG_n' : 1.23,
        'vG_m' : (1-1/1.23),
        'vG_a' : 50000,
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
        }
    
        
class Soil_domain:
    def __init__ (self):
        self.Depth_top = []
        self.Depth_bottom = []
        self.precision = []
        self.node_initial = []




