# -------------------------------------------
# General Input
# -------------------------------------------
site = "small_strain"
tstart = "2024-01-01"
tend = "2025-01-01"


# -------------------------------------------
# Groundwater
# -------------------------------------------
gwli = 0 		#initial groundwater level
Ref_Time = "10:00" 	#Reference Time
GL_to_NAP = -1.48	#[m]

# -------------------------------------------
# Soil depths and types
# -------------------------------------------
depths = [-1, -11.0, -12.0]
types = ['aquifer', 'confining_layer', 'aquifer']

# -------------------------------------------
# Boundary condition type at Top
# -------------------------------------------
#0 is Dirichlet boundary condition. (Pressure [kPa])
#1 is Neumann boundary condition (flux [m/day])
tp_bc = 0 #Top boundary condition type

bc = 0 # if bottom boudary condition is given as time series bc = 1