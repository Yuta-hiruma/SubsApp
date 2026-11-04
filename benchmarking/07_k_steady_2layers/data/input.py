# -------------------------------------------
# General Input
# -------------------------------------------
site = "benchmark"
tstart = "2024-01-01"
#tend = "2024-10-10"
tend = "2025-05-15"
#tend = "2026-09-27"


# -------------------------------------------
# Groundwater
# -------------------------------------------
gwli = -50 		#initial groundwater level
Ref_Time = "10:00" 	#Reference Time
GL_to_NAP = -10.0	#[m]

# -------------------------------------------
# Soil depths and types
# -------------------------------------------
depths = [-0.5,-4.0]
types = ['Clay1','Clay2']


# -------------------------------------------
# Boundary condition type at Top
# -------------------------------------------
#0 is Dirichlet boundary condition. (Pressure [kPa])
#1 is Neumann boundary condition (flux [m/day])
tp_bc = 1 #Top boundary condition type

bc = 0 # if bottom boudary condition is given as time series bc = 1