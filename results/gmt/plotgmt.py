import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib.colors import LinearSegmentedColormap



gmt_dir = "gmt"
#_/_/_/_/_/_/_/_/_/_/_/_/
#Se
#_/_/_/_/_/_/_/_/_/_/_/_/
# Read and reshape the data
Se_vG_file = "../Se_vG.dat"
Se_vG = pd.read_csv(
    Se_vG_file,
    sep=r"\s+",
    skiprows=1,
    names=["time", "depth", "effective_saturation"],
)

Se_grid = (
    Se_vG.pivot(
        index="depth",
        columns="time",
        values="effective_saturation",
    )
    .sort_index()
    .sort_index(axis=1)
)

simulation_time = Se_grid.columns.to_numpy(dtype=float)
depth = Se_grid.index.to_numpy(dtype=float)
Se_values = Se_grid.to_numpy(dtype=float)

# GMT xyz2grd requires regularly spaced coordinates
if len(simulation_time) < 2 or len(depth) < 2:
    raise ValueError("At least two time and depth coordinates are required.")

dx = np.diff(simulation_time)
dy = np.diff(depth)

if not np.allclose(dx, dx[0]):
    raise ValueError("Time coordinates must be regularly spaced for xyz2grd.")

if not np.allclose(dy, dy[0]):
    raise ValueError("Depth coordinates must be regularly spaced for xyz2grd.")

# Export time, depth, effective saturation as GMT x/y/z data
time_mesh, depth_mesh = np.meshgrid(simulation_time, depth)
xyz = np.column_stack((time_mesh.ravel(), depth_mesh.ravel(), Se_values.ravel()))

np.savetxt("Se_vG_xyz.dat", xyz, fmt="%.12g")

# Preserve the original range and threshold for the groundwater boundary
Se_min = np.nanmin(Se_values)
Se_max = np.nanmax(Se_values)

if Se_min == Se_max:
    raise ValueError(
        "Effective saturation is constant; a non-zero color range is required."
    )

bounds = np.linspace(Se_min, Se_max, 31)
if Se_min >= 1.0:
    bounds = np.linspace(0.0, max(1.0, Se_max), 31)

boundary_Se = bounds[-2]

# Create a smooth GMT CPT using 256 interpolated color segments
n_cpt = 256
cpt_bounds = np.linspace(bounds[0], bounds[-1], n_cpt + 1)

colors = plt.get_cmap("Blues")(np.linspace(0, 1, n_cpt + 1))
colors[0] = (1, 1, 1, 1)  # White at the lowest CPT value

def rgb_string(rgba):
    rgb = np.round(np.asarray(rgba[:3]) * 255).astype(int)
    return f"{rgb[0]}/{rgb[1]}/{rgb[2]}"

cpt_file = "Se_vG.cpt"
with open(cpt_file, "w", encoding="utf-8") as cpt:
    for i in range(n_cpt):
        cpt.write(
            f"{cpt_bounds[i]:.12g} {rgb_string(colors[i])} "
            f"{cpt_bounds[i + 1]:.12g} {rgb_string(colors[i + 1])}\n"
        )
    cpt.write(f"B {rgb_string(colors[0])}\n")
    cpt.write(f"F {rgb_string(colors[-1])}\n")
    cpt.write("N 128/128/128\n")

# Write the outer grid bounds and grid spacing
west = simulation_time[0] - dx[0] / 2
east = simulation_time[-1] + dx[0] / 2
south = depth[0] - dy[0] / 2
north = depth[-1] + dy[0] / 2

Path("Se_vG_gridinfo.txt").write_text(
    f"{west:.12g} {east:.12g} {south:.12g} {north:.12g} "
    f"{dx[0]:.12g} {dy[0]:.12g}\n",
    encoding="utf-8",
)

# Calculate and export the simulated groundwater-level boundary
boundary_depth = np.full(len(simulation_time), np.nan)

for j in range(len(simulation_time)):
    profile = Se_values[:, j]

    crossings = np.where(
        ((profile[:-1] < boundary_Se) & (profile[1:] >= boundary_Se))
        | ((profile[:-1] >= boundary_Se) & (profile[1:] < boundary_Se))
    )[0]

    if crossings.size:
        crossing_depths = (
            depth[crossings]
            + (boundary_Se - profile[crossings])
            / (profile[crossings + 1] - profile[crossings])
            * (depth[crossings + 1] - depth[crossings])
        )
        boundary_depth[j] = crossing_depths[
            np.argmin(np.abs(crossing_depths))
        ]

np.savetxt(
    "groundwater_level.dat",
    np.column_stack((simulation_time, boundary_depth)),
    fmt="%.12g",
)
#_/_/_/_/_/_/_/_/_/_/_/_/
#theta
#_/_/_/_/_/_/_/_/_/_/_/_/

# Read and reshape the data
theta_file = "../theta.dat"
theta = pd.read_csv(
    theta_file,
    sep=r"\s+",
    skiprows=1,
    names=["time", "depth", "theta"],
)

theta_grid = (
    theta.pivot(
        index="depth",
        columns="time",
        values="theta",
    )
    .sort_index()
    .sort_index(axis=1)
)

simulation_time = theta_grid.columns.to_numpy(dtype=float)
depth = theta_grid.index.to_numpy(dtype=float)
theta_values = theta_grid.to_numpy(dtype=float)

# GMT xyz2grd requires regularly spaced coordinates
if len(simulation_time) < 2 or len(depth) < 2:
    raise ValueError("At least two time and depth coordinates are required.")

dx = np.diff(simulation_time)
dy = np.diff(depth)

if not np.allclose(dx, dx[0]):
    raise ValueError("Time coordinates must be regularly spaced for xyz2grd.")

if not np.allclose(dy, dy[0]):
    raise ValueError("Depth coordinates must be regularly spaced for xyz2grd.")

# Export time, depth, effective saturation as GMT x/y/z data
time_mesh, depth_mesh = np.meshgrid(simulation_time, depth)
xyz = np.column_stack((time_mesh.ravel(), depth_mesh.ravel(), theta_values.ravel()))

np.savetxt("theta_xyz.dat", xyz, fmt="%.12g")

# Preserve the original range and threshold for the groundwater boundary
theta_min = np.nanmin(theta_values)
theta_max = np.nanmax(theta_values)

if theta_min == theta_max:
    raise ValueError(
        "Effective saturation is constant; a non-zero color range is required."
    )

bounds = np.linspace(theta_min, theta_max, 31)


boundary_theta = bounds[-2]

# Create a smooth GMT CPT using 256 interpolated color segments
n_cpt = 256
cpt_bounds = np.linspace(bounds[0], bounds[-1], n_cpt + 1)

colors = plt.get_cmap("YlGnBu")(np.linspace(0, 1, n_cpt + 1))
colors[0] = (1, 1, 1, 1)  # White at the lowest CPT value

def rgb_string(rgba):
    rgb = np.round(np.asarray(rgba[:3]) * 255).astype(int)
    return f"{rgb[0]}/{rgb[1]}/{rgb[2]}"

cpt_file = "theta.cpt"
with open(cpt_file, "w", encoding="utf-8") as cpt:
    for i in range(n_cpt):
        cpt.write(
            f"{cpt_bounds[i]:.12g} {rgb_string(colors[i])} "
            f"{cpt_bounds[i + 1]:.12g} {rgb_string(colors[i + 1])}\n"
        )
    cpt.write(f"B {rgb_string(colors[0])}\n")
    cpt.write(f"F {rgb_string(colors[-1])}\n")
    cpt.write("N 128/128/128\n")

# Write the outer grid bounds and grid spacing
west = simulation_time[0] - dx[0] / 2
east = simulation_time[-1] + dx[0] / 2
south = depth[0] - dy[0] / 2
north = depth[-1] + dy[0] / 2

Path("theta_gridinfo.txt").write_text(
    f"{west:.12g} {east:.12g} {south:.12g} {north:.12g} "
    f"{dx[0]:.12g} {dy[0]:.12g}\n",
    encoding="utf-8",
)


#_/_/_/_/_/_/_/_/_/_/_/_/
#pressure
#_/_/_/_/_/_/_/_/_/_/_/_/
# Read and reshape the data
pressure_file = Path("../Pressure.dat")
pressure = pd.read_csv(
    pressure_file,
    sep=r"\s+",
    skiprows=1,
    names=["time", "depth", "pressure"],
)

pressure_grid = (
    pressure.pivot(
        index="depth",
        columns="time",
        values="pressure",
    )
    .sort_index()
    .sort_index(axis=1)
)

simulation_time = pressure_grid.columns.to_numpy(dtype=float)
depth = pressure_grid.index.to_numpy(dtype=float)
pressure_values = pressure_grid.to_numpy(dtype=float)

# GMT xyz2grd requires regularly spaced coordinates
if len(simulation_time) < 2 or len(depth) < 2:
    raise ValueError("At least two time and depth coordinates are required.")

dx = np.diff(simulation_time)
dy = np.diff(depth)

if not np.allclose(dx, dx[0]):
    raise ValueError("Time coordinates must be regularly spaced for xyz2grd.")

if not np.allclose(dy, dy[0]):
    raise ValueError("Depth coordinates must be regularly spaced for xyz2grd.")

# Export time, depth, and pressure as GMT x/y/z data
time_mesh, depth_mesh = np.meshgrid(simulation_time, depth)
xyz = np.column_stack(
    (time_mesh.ravel(), depth_mesh.ravel(), pressure_values.ravel())
)
np.savetxt("pressure_xyz.dat", xyz, fmt="%.12g")

pressure_min = np.nanmin(pressure_values)
pressure_max = np.nanmax(pressure_values)

cpt_min = min(pressure_min, 0.0)
cpt_max = max(pressure_max, 0.0)

if cpt_min == cpt_max:
    raise ValueError("All pressure values are zero; a non-zero CPT range is required.")

# Select the relevant part of the diverging scale
if cpt_min < 0 < cpt_max:
    zero_position = (0.0 - cpt_min) / (cpt_max - cpt_min)
    color_stops = [
        (0.0, "red"),
        (zero_position, "white"),
        (1.0, "blue"),
    ]
elif cpt_min == 0:
    # Non-negative range: white at zero, blue toward the maximum
    color_stops = [(0.0, "white"), (1.0, "blue")]
else:
    # Non-positive range: red at the minimum, white at zero
    color_stops = [(0.0, "red"), (1.0, "white")]

cmap = LinearSegmentedColormap.from_list(
    "red_white_blue",
    color_stops,
    N=257,
)

# Include zero as an exact CPT boundary
n_cpt = 256
cpt_bounds = np.unique(
    np.append(np.linspace(cpt_min, cpt_max, n_cpt + 1), 0.0)
)
colors = cmap((cpt_bounds - cpt_min) / (cpt_max - cpt_min))

def rgb_string(rgba):
    rgb = np.round(np.asarray(rgba[:3]) * 255).astype(int)
    return f"{rgb[0]}/{rgb[1]}/{rgb[2]}"

# Write a continuous GMT CPT
with open("pressure.cpt", "w", encoding="utf-8") as cpt:
    for i in range(len(cpt_bounds) - 1):
        cpt.write(
            f"{cpt_bounds[i]:.12g} {rgb_string(colors[i])} "
            f"{cpt_bounds[i + 1]:.12g} {rgb_string(colors[i + 1])}\n"
        )
    cpt.write(f"B {rgb_string(colors[0])}\n")
    cpt.write(f"F {rgb_string(colors[-1])}\n")
    cpt.write("N 128/128/128\n")

# Write the outer grid bounds and grid spacing
west = simulation_time[0] - dx[0] / 2
east = simulation_time[-1] + dx[0] / 2
south = depth[0] - dy[0] / 2
north = depth[-1] + dy[0] / 2

Path("pressure_gridinfo.txt").write_text(
    f"{west:.12g} {east:.12g} {south:.12g} {north:.12g} "
    f"{dx[0]:.12g} {dy[0]:.12g}\n",
    encoding="utf-8",
)