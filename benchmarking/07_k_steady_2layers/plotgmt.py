import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

CASE_DIR = Path(__file__).resolve().parent
BENCH_OUTPUT_DIR = CASE_DIR / "output"
RESULTS_DIR = CASE_DIR / "output"

output_dir = BENCH_OUTPUT_DIR
gmt_dir = output_dir / "gmt"
gmt_dir.mkdir(parents=True, exist_ok=True)

# Read and reshape the data
Se_vG_file = RESULTS_DIR / "Se_vG.dat"
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

np.savetxt(gmt_dir / "Se_vG_xyz.dat", xyz, fmt="%.12g")

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

cpt_file = gmt_dir / "Se_vG.cpt"
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

(gmt_dir / "Se_vG_gridinfo.txt").write_text(
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
    gmt_dir / "groundwater_level.dat",
    np.column_stack((simulation_time, boundary_depth)),
    fmt="%.12g",
)