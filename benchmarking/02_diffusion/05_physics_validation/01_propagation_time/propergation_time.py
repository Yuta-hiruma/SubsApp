import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.interpolate import interp1d
import os

# Parameters
z_f = -1.35
z_f_2 = - 0.4  # depth of profile [m]
theta_s = 0.53
theta_r = 0.0
theta_top = 0.38
#theta =[0.381, 0.385,0.388,0.392,0.396,0.400,0.404,0.408,0.413,0.418,0.423,0.428,0.433,0.439,0.445,0.452,0.458,0.465,0.473,0.481,0.490,0.499,0.508,0.518,0.528
#,0.530,0.530]
n = 1.28
m = 1 - 1 / n
_lambda = 0.5
#K0_list = [1,0.5,0.1,0.02,0.01] 
K0_list = [0.01,0.05,0.1,0.5,1]# m/day
print("Ko=", K0_list)  
cell = 50
#cell = len(theta)  #cell number
print(cell)


# Depth range
z_vals = np.linspace(0, abs(z_f), cell)
z_mids = 0.5 * (z_vals[:-1] + z_vals[1:])

cell_gaps = np.diff(z_vals)
dz= cell_gaps[1]
# Containers
Se_profiles = {}
K_profiles = {}
T_profiles = {}


def theta_z(z):
    return (theta_s - theta_top) / abs(z_f) * z + theta_top
    #return (theta_s - theta_top) / abs(z_f_2) * z + theta_top

def Se_z(z):
    Se = (theta_z(z) - theta_r) / (theta_s - theta_r)
    if Se >=1:
        Se = 1
    return Se

def K_z(z, K0):
    Se = Se_z(z)
    #if Se <= 0:
     #   return 1e-12
    return (K0 * Se**_lambda * (1 - (1 - Se**(1 / m))**m)**2)/dz

def integrand(z, K0):
    K = K_z(z, K0)
    return  K



# Compute for each K0
for K0 in K0_list:
    Se_vals = [Se_z(z) for z in z_mids]
    #print(Se_vals)
    K_vals = [K_z(z, K0) for z in z_mids]
    T_vals = []
    integral = 0
    for z in z_mids:

        integral = integral + dz/integrand(z, K0)
        T_vals.append(integral)

    Se_profiles[K0] = Se_vals
    K_profiles[K0] = K_vals
    T_profiles[K0] = T_vals



#Out put----------------------------------------------
output_dir = "output"
os.makedirs(output_dir, exist_ok=True)

# Save each profile to a separate .dat file
for K0 in K0_list:
    data = np.column_stack((-z_mids, Se_profiles[K0], K_profiles[K0], T_profiles[K0]))
    header = "Depth[m]  Se[-]  K[m/day]  T[day]"
    filename = os.path.join(output_dir, f"profile_K0_{K0}.dat")
    np.savetxt(filename, data, header=header, fmt="%.5f", delimiter="\t")


#Plot ------------------------------------------------
# Se vs z ---
plt.figure(figsize=(6, 4))

plt.plot(Se_profiles[K0], -z_mids)
plt.xlabel('Effective Saturation $S_e(z)$ [-]')
plt.ylabel('Depth $z$ [m]')
plt.title('Effective Saturation vs Depth')
plt.grid(True)
plt.tight_layout()
plt.savefig('Se.png')


# K vs z ---
plt.figure(figsize=(6, 4))
for K0 in K0_list:
    plt.plot(K_profiles[K0], -z_mids, label=f'$K_0$={K0} m/day')
plt.xlabel('Hydraulic Conductivity $K(z)$  * unit lenght $dz$ [day]')
plt.ylabel('Depth $z$ [m]')
plt.xticks(rotation=45)
plt.subplots_adjust(left=0.2, bottom=0.2) 
plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
plt.title('Hydraulic Conductivity vs Depth')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('K_multiple.png')


# T(z) vs z ---
plt.figure(figsize=(6, 4))
for K0 in K0_list:
    plt.plot(T_profiles[K0], -z_mids, label=f'$K_0$={K0} m/day')
plt.xlabel('Diffusion Time $T(z)$ [day]')
plt.ylabel('Depth $z$ [m]')
plt.xticks(rotation=45)
plt.subplots_adjust(left=0.2, bottom=0.2) 
plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
plt.title('Diffusion Time (ΔP = 10kPa)')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('Tz_multiple.png')

# K vs Se ---
plt.figure(figsize=(6, 4))
for K0 in K0_list:
    plt.plot(Se_profiles[K0], K_profiles[K0], label=f'$K_0$={K0} m/day')

plt.xlabel('Effective Saturation $S_e$ [-]')
plt.ylabel('Hydraulic Conductivity $K$ [m/day]')
plt.title('Hydraulic Conductivity vs Effective Saturation')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('K_vs_Se.png', dpi=300)
