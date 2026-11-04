# -------------------------------------------
# General Input
# -------------------------------------------
site = "nieuwolda"
tstart = "2022-07-01"
tend = "2024-07-31"



# -------------------------------------------
# Groundwater
# -------------------------------------------
gwli = -156.5 		#initial groundwater level
Ref_Time = "10:00" 	#Reference Time
GL_to_NAP = -1.48	#[m]

# -------------------------------------------
# Soil depths and types
# -------------------------------------------
depths = [-2.7, -5.0, -5.70,-6.0, -6.50]
types = ["Clay_top", "Peat1", "Clay2", "Peat2" , "Sand"] 


# -------------------------------------------
# Pressure boundary condition type
# -------------------------------------------
#0 is Dirichlet boundary condition. (Pressure [kPa])
#1 is Neumann boundary condition (flux [m/day])
tp_bc = 0
bt_bc = 0