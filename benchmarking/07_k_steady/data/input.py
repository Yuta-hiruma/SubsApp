# -------------------------------------------
# General Input
# -------------------------------------------
site = "benchmark"
tstart = "2024-01-01"
tend = "2024-04-10"


# -------------------------------------------
# Groundwater
# -------------------------------------------
gwli = -50 		#initial groundwater level
Ref_Time = "10:00" 	#Reference Time
GL_to_NAP = -10.0	#[m]

# -------------------------------------------
# Soil depths and types
# -------------------------------------------
depths = [-1.0]
types = ['Clay']


# -------------------------------------------
# Pressure boundary condition type
# -------------------------------------------
#0 is Dirichlet boundary condition. (Pressure [kPa])
#1 is Neumann boundary condition (flux [m/day])
tp_bc = 0
bt_bc = 0