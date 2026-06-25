

Soil_parameters ={
    'Clay':{
        'Material_parameters':{
         'gamma_s' : 2.31,
        'gamma_w' : 1,
        },
        
        'Compression': {
        'kappa': 0.078,
        'lambda': 0.114,
        'mu': 0.045,
        'Tau_ref': 1,
        'Pg': 400,
        'OCR': 1.2,
        'c' : 10,
        'phi' : 20
        },
        
        'Water_retention_curve':{
        'vG_WCS' : 0.89,
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
        'K_tot' : 0.1,
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

