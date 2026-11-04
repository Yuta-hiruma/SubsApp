# -------------------------------------------
# General Input
# -------------------------------------------
site = "benchmark"
tstart = "2024-01-01"
tend = "2024-01-10"



# -------------------------------------------
# Groundwater
# -------------------------------------------
gwli = -50 		#initial groundwater level
Ref_Time = "10:00" 	#Reference Time
GL_to_NAP = -0.05	#[m]

# -------------------------------------------
# Soil depths and types
# -------------------------------------------
depths = [-1.0]
types = ['Clay']


# -------------------------------------------
# Boundary condition type at Top
# -------------------------------------------
#0 is Dirichlet boundary condition. (Pressure [kPa])
#1 is Neumann boundary condition (flux [m/day])
tp_bc = 0 #Top boundary condition type

bc = 0 # if bottom boudary condition is given as time series bc = 1