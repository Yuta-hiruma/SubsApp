import numpy as np
import pandas as pd
import pygmt

#Reading input data ----------
import sys
import os
current_dir = os.getcwd()
parent_dir = os.path.abspath("../data")
sys.path.insert(0, parent_dir)
from Input_soil_material_parameters import *
from input import *
sys.path.pop(0)
os.chdir(current_dir)

k = pd.read_csv(
    "k.dat", sep=r"\s+", skiprows=1,
    names=["time", "depth", "K"]
)
Se = pd.read_csv(
    "Se_vG.dat", sep=r"\s+", skiprows=1,
    names=["time", "depth", "effective_saturation"]
)

nodes = [-0.275, -0.325, -0.375]
#times = [8, 9]
times = [9]


K_sat = Soil_parameters['Clay']["Permeability_curve_Perm"]["K_tot"]
tau = Soil_parameters['Clay']["Permeability_curve_Perm"]["tau"]
m = Soil_parameters['Clay']["Water_retention_curve"]["vG_m"]


def get_value(data, time, depth, column):
    """Return the value matching a time and depth, or None if not found."""
    rows = data[
        (data["time"] == time)
        & np.isclose(data["depth"], depth)
    ]
    return None if rows.empty else rows[column].iloc[0]


def update_K_Kmin (S_e,K_sat,m,tau): 
    K_out = K_sat*(S_e**tau)*(1-(1-S_e**(1/m)**m))**2
    return (K_out)  


results = []

for time in times:
    # Print K and saturation at each selected node.
    for depth in nodes:
        k_value = get_value(k, time, depth, "K")
        se_value = get_value(Se, time, depth, "effective_saturation")

        if k_value is None or se_value is None:
            print(f"No data found for time={time}, depth={depth}")
        else:
            print(f"time={time}, depth={depth}: K={k_value}, Se={se_value}")

    # Calculate means between each pair of adjacent nodes.
    for depth_1, depth_2 in zip(nodes[:-1], nodes[1:]):
        k_1 = get_value(k, time, depth_1, "K")
        k_2 = get_value(k, time, depth_2, "K")
        Se_1 = get_value(Se, time, depth_1, "effective_saturation")
        Se_2 = get_value(Se, time, depth_2, "effective_saturation")
        print("Se_1", Se_1)
        print("Se_1", Se_2)
        if k_1 is None or k_2 is None:
            print(f"Cannot calculate midpoint K for time={time}, "
                  f"depths={depth_1}, {depth_2}")
            continue

        if k_1 <= 0 or k_2 <= 0:
            raise ValueError("K values must be positive to calculate the logarithmic mean.")

        midpoint = (depth_1 + depth_2) / 2

        Se_ave = (Se_2+Se_1)/2
        print("Se_ave",Se_ave)
        k_se = update_K_Kmin(Se_ave,K_sat,m,tau) 

        # Harmonic mean
        k_harmonic = 2 * k_1 * k_2 / (k_1 + k_2)

        # Logarithmic mean; use its limit when the two values are equal.
        if np.isclose(k_1, k_2):
            k_logarithmic = k_1
        else:
            k_logarithmic = (k_2 - k_1) / np.log(abs(k_2 / k_1))

        results.append({
            "time": time,
            "depth_midpoint": midpoint,
            "K_harmonic": k_harmonic,
            "K_logarithmic": k_logarithmic,
            "K_se": k_se,
        })

results = pd.DataFrame(results)

fig = pygmt.Figure()

colors = ["royalblue", "firebrick"]
depth_pad = 0.02

# Set the plot region to include the original and calculated K values.
x_values = []

for time in times:
    original = k[k["time"] == time]
    original = original[
        original["depth"].apply(
            lambda depth: any(np.isclose(depth, node) for node in nodes)
        )
    ]
    subset = results[results["time"] == time]

    x_values.extend(original["K"].dropna().tolist())
    for column in ["K_harmonic", "K_logarithmic", "K_se"]:
        x_values.extend(subset[column].dropna().tolist())

xmin, xmax = min(x_values), max(x_values)
x_pad = 0.05 * (xmax - xmin) if xmax > xmin else 0.05 * abs(xmax)

fig.basemap(
    region=[
        xmin - x_pad,
        xmax + x_pad,
        min(nodes) - depth_pad,
        max(nodes) + depth_pad,
    ],
    projection="X12c/9c",
    frame=[
        "xaf+lHydraulic conductivity, K [m/day]",
        "yaf+lDepth [m]",
        "WSen",
    ],
)

def add_value_labels(values, depths, prefix, color, offset, justify):
    fig.text(
        x=values,
        y=depths,
        text=[f"{prefix}={value:.6f}" for value in values],
        font=f"7p,Helvetica,{color}",
        justify=justify,
        offset=offset,
    )


for i, time in enumerate(times):
    color = colors[i]

    # Original K values at the nodes: connected line and black circles.
    original = k[k["time"] == time].copy()
    original = original[
        original["depth"].apply(
            lambda depth: any(np.isclose(depth, node) for node in nodes)
        )
    ].sort_values("depth")

    fig.plot(
        x=original["K"],
        y=original["depth"],
        pen="1p,black",
    )
    fig.plot(
        x=original["K"],
        y=original["depth"],
        style="c0.12c",
        fill="black",
        pen="0.5p,black",
        label=f"K at nodes, t = {time} day",
    )

    # Calculated values at midpoint depths: symbols without connecting lines.
    subset = results[results["time"] == time]

    fig.plot(
        x=subset["K_harmonic"],
        y=subset["depth_midpoint"],
        style="s0.15c",
        fill=color,
        pen=f"0.5p,{color}",
        label=f"Harmonic mean, t = {time} day",
    )
    add_value_labels(
        subset["K_harmonic"].to_numpy(),
        subset["depth_midpoint"].to_numpy(),
        "H", color, "0.12c/0.10c", "BL"
    )

    fig.plot(
        x=subset["K_logarithmic"],
        y=subset["depth_midpoint"],
        style="t0.15c",
        fill=color,
        pen=f"0.5p,{color}",
        label=f"Logarithmic mean, t = {time} day",
    )
    add_value_labels(
        subset["K_logarithmic"].to_numpy(),
        subset["depth_midpoint"].to_numpy(),
        "L", color, "0.12c/0c", "ML"
    )

    fig.plot(
        x=subset["K_se"],
        y=subset["depth_midpoint"],
        style="d0.15c",
        fill=color,
        pen=f"0.5p,{color}",
        label=f"K(Se), t = {time} day",
    )
    add_value_labels(
        subset["K_se"].to_numpy(),
        subset["depth_midpoint"].to_numpy(),
        "K(Se)", color, "0.12c/-0.10c", "TL"
    )

fig.legend()
fig.savefig("k_interval.png")