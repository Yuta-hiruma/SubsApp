# -------------------------------------------
# General Input
# -------------------------------------------
site = "diffusion"
tstart = "2024-01-01"
tend = "2024-08-20"



# -------------------------------------------
# Groundwater
# -------------------------------------------
gwli = -156.5 		#initial groundwater level
Ref_Time = "10:00" 	#Reference Time
GL_to_NAP = -1.48	#[m]

# -------------------------------------------
# Soil depths and types
# -------------------------------------------
depths = [-6.0, -6.50]
types = ['Clay_top_SWP',  'Sand_SWP']

# -------------------------------------------
# Boundary condition type at Top
# -------------------------------------------
#0 is Dirichlet boundary condition. (Pressure [kPa])
#1 is Neumann boundary condition (flux [m/day])
tp_bc = 0 #Top boundary condition type

bc = 0 # if bottom boudary condition is given as time series bc = 1