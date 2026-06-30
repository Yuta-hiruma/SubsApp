

Soil_parameters ={
    'Peat1_SWP':{
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
        'OCR': 1.5,
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

Soil_parameters['Peat2_SWP'] = {
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
        'OCR': 1.5,
        'c' : 10,
        'phi' : 20
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
        'K_tot' : 0.1,
        'K_min' : 0.0001,
        'tau' : 0.5
        }      
        
        }


 
Soil_parameters ['Clay_top_SWP']= {      
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
        'OCR': 1.5,
        'c' : 10,
        'phi' : 20
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
        'K_tot' : 0.01,
        'K_min' : 0.0001,
        'tau' : 0.5
        }
        
        }
Soil_parameters ['Clay2_SWP']= {      
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
        'OCR': 1.5,
        'c' : 10,
        'phi' : 20
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
        'K_tot' : 0.01,
        'K_min' : 0.0001,
        'tau' : 0.5
        }
        
        }
Soil_parameters ['Sand_SWP']= {      
        'Material_parameters':{
        'gamma_s' : 2.31,
        'gamma_w' : 1,
        },
        
        'Compression':{
        'kappa' : 0.0868,   
        'lambda'  : 0.1302,
        'mu'  : 0.0028229,
        'Tau_ref'  : 1,
        'Pg' : 400,
        'OCR' : 10,
        'c' : 10,
        'phi' : 20
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
        'K_tot' : 0.01,
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

