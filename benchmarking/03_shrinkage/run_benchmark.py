from pathlib import Path
import shutil
import subprocess
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib import ticker
#Reading input data ----------
import sys
import os
current_dir = os.getcwd()
parent_dir = os.path.abspath("data")
sys.path.insert(0, parent_dir)
from Input_soil_material_parameters import *
sys.path.pop(0)
os.chdir(current_dir)


CASE_DIR = Path(__file__).resolve().parent
ROOT_DIR = CASE_DIR.parents[1]

SRC_DIR = ROOT_DIR / "src"
DATA_DIR = ROOT_DIR / "data"
RESULTS_DIR = ROOT_DIR / "results"
GMT_DIR = ROOT_DIR / "results/gmt"

BENCH_DATA_DIR = CASE_DIR / "data"
BENCH_INPUT_DIR = CASE_DIR / "input"
BENCH_OUTPUT_DIR = CASE_DIR / "output"

BACKUP_DIR = CASE_DIR / "_backup_runtime_files"


#Read measurement data
measure_void_moist = "measurement/void_moisture.dat"
measure_suction_time = "measurement/suction_time.dat"
void_moist_tp = pd.read_csv(measure_void_moist, sep=r"\s+", skiprows=1, names=["moist","void"])
suction_time_bp = pd.read_csv(measure_suction_time, sep=r"\s+", skiprows=1, names=["time", "suction"])


def backup_path(path: Path):
    backup = BACKUP_DIR / path.relative_to(ROOT_DIR)
    backup.parent.mkdir(parents=True, exist_ok=True)

    if path.exists():
        if path.is_dir():
            shutil.copytree(path, backup, dirs_exist_ok=True)
        else:
            shutil.copy2(path, backup)


def restore_path(path: Path):
    backup = BACKUP_DIR / path.relative_to(ROOT_DIR)

    if path.exists():
        if path.is_dir():
            shutil.rmtree(path)
        else:
            path.unlink()

    if backup.exists():
        if backup.is_dir():
            shutil.copytree(backup, path, dirs_exist_ok=True)
        else:
            path.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(backup, path)


def copy_case_files():
    DATA_DIR.mkdir(exist_ok=True)

    shutil.copy2(BENCH_DATA_DIR / "input.py", DATA_DIR / "input.py")
    shutil.copy2(
        BENCH_DATA_DIR / "Input_soil_material_parameters.py",
        DATA_DIR / "Input_soil_material_parameters.py",
    )

    if (SRC_DIR / "input").exists():
        shutil.rmtree(SRC_DIR / "input")

    shutil.copytree(BENCH_INPUT_DIR, SRC_DIR / "input")


def run_model():
    RESULTS_DIR.mkdir(exist_ok=True)

    subprocess.run(
        [sys.executable, "main.py"],
        cwd=SRC_DIR,
        check=True,
    )


def collect_results():
    BENCH_OUTPUT_DIR.mkdir(exist_ok=True)

    Pressure_fig = RESULTS_DIR /  "gmt/pressure.png"
    Se_vG_fig = RESULTS_DIR /  "gmt/Se_vG.png"

    shutil.copy2(Pressure_fig, BENCH_OUTPUT_DIR / "pressure.png",)
    shutil.copy2(Se_vG_fig, BENCH_OUTPUT_DIR / "Se_vG.png",)



def plotting():
    # read measurement data-----
    measure_void_moist = "measurement/void_moisture.dat"
    measure_suction_time = "measurement/suction_time.dat"
    void_moist_tp = pd.read_csv(measure_void_moist, sep=r"\s+", skiprows=1, names=["moist","void"])
    suction_time_bp = pd.read_csv(measure_suction_time, sep=r"\s+", skiprows=1, names=["time", "suction"])
    
    # read simulation result-----
    stress_strain_file = RESULTS_DIR / "stress_strain.dat"
    pressure_void_file = RESULTS_DIR / "pressure_void_simulation.dat"
    void_time_file = RESULTS_DIR / "void_time_simulation.dat"
    pressure_time_file = RESULTS_DIR / "pressure_time_simulation.dat"
    sscc_MR_file = RESULTS_DIR / "sscc_MR_simulation.dat"
    #Se_vG_file = RESULTS_DIR / "Se_vG.dat"
    stress_strain = pd.read_csv(stress_strain_file,sep=r"\s+",skiprows=1,names=["sigma_eff", "strain"])
    pressure_void = pd.read_csv(pressure_void_file,sep=r"\s+",skiprows=1,names=["pressure", "void"])
    void_time = pd.read_csv(void_time_file,sep=r"\s+",skiprows=1,names=["time", "void"])
    pressure_time = pd.read_csv(pressure_time_file,sep=r"\s+",skiprows=1,names=["time", "pressure"])
    sscc_MR = pd.read_csv(sscc_MR_file,sep=r"\s+",skiprows=1,names=["moisture_ratio", "void"])
    #Se_vG = pd.read_csv(Se_vG_file,sep=r"\s+",skiprows=1,names=["time", "depth", "effective_saturation"],)
    #Se_grid = Se_vG.pivot(index="depth",columns="time",values="effective_saturation",).sort_index().sort_index(axis=1)
    

    #-----------------------------------
    #Effective stress - strain
    plt.figure()
    plt.plot(stress_strain["sigma_eff"],stress_strain["strain"],marker="s",linestyle=":",color="black")
    plt.title("Stress strain")
    plt.xlabel(r"Effective stress $\sigma'$ [kPa]")
    plt.ylabel(r"Strain $\varepsilon$ [-]")
    plt.xscale("log")
    plt.grid(True,which="both",linestyle="--",alpha=0.5)
    plt.tight_layout()
    textstr = (
    rf"$\kappa = {Soil_parameters["Clay"]["Compression"]["kappa"]:.3f}$" "\n"
    rf"$\lambda = {Soil_parameters["Clay"]["Compression"]["lambda"]:.3f}$" "\n"
    rf"$\mu = {Soil_parameters["Clay"]["Compression"]["mu"]:.3f}$"
    )
    #plt.text(0.05, 0.05,textstr,transform=plt.gca().transAxes,fontsize=11,verticalalignment="bottom",
    #bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))
    plt.savefig(BENCH_OUTPUT_DIR / "stress_strain.png",dpi=300)
    plt.close()

    #-----------------------------------
    #Pore pressure - void ratio
    plt.figure()
    plt.plot(pressure_void["pressure"],pressure_void["void"],marker="s",linestyle=":",label="Simulation")
    plt.title("Pressure - void ratio")
    plt.xlabel(r"Pore pressure $\psi$ [kPa]")
    plt.ylabel(r"Void ratio $e$ [-]")
    plt.legend()
    plt.grid(True,which="both",linestyle="--",alpha=0.5)
    plt.tight_layout()
    plt.savefig(BENCH_OUTPUT_DIR / "pressure_void.png",dpi=300)
    plt.close()

    #-----------------------------------
    #Void ratio - time
    plt.figure()
    #plt.plot(suction_time_bp["time"],void_suction_bp["void"],".",color="black",label="Measurement")
    plt.plot(void_time["time"],void_time["void"],label="Simulation")
    plt.title("Void ratio - time")
    plt.xlabel("Time [day]")
    plt.ylabel(r"Void ratio $e$")
    plt.legend()
    plt.grid(True,which="both",linestyle="--",alpha=0.5)
    plt.tight_layout()
    plt.savefig(BENCH_OUTPUT_DIR / "void_time.png",dpi=300)
    plt.close()
    #-----------------------------------
    #Soil shrinkage characteristic curve
    plt.figure()
    plt.plot(void_moist_tp["moist"],void_moist_tp["void"],".",markeredgewidth=0.1,label="Measurement")
    plt.plot(sscc_MR["moisture_ratio"],sscc_MR["void"],marker="s",linestyle=":",color="black",label="Simulation by SubsApp")
    x = void_moist_tp["void"]
    plt.plot(x,x,color="black", ls="--" ,label="Saturation")
    # Position the equation below the line
    x_pos = x.min() + 0.55 * (x.max() - x.min())
    offset = 0.05 * (x.max() - x.min())
    plt.text(x_pos,x_pos - offset,r"$e = \vartheta = \frac{V_v}{V_s}, \quad V_v = V_w$",rotation=45,rotation_mode="anchor",ha="center",va="top",color="black")
    #plt.title("Soil shrinkage characteristic curve")
    plt.xlabel(r"Moisture ratio $\vartheta$")
    plt.ylabel(r"Void ratio $e$")
    plt.legend()
    plt.grid(True,which="both",linestyle="--",alpha=0.5)

    plt.tight_layout()
    plt.savefig(BENCH_OUTPUT_DIR / "sscc_MR.png",dpi=300)
    plt.close()

def run_gmt_plots():
    commands = [
        ["python3", "plotgmt.py"],
        ["bash", "plot_Se_vG.sh"],
        ["bash", "plot_theta.sh"],
        ["bash", "plot_pressure.sh"],
    ]

    for command in commands:
        subprocess.run(command, cwd=GMT_DIR, check=True)
    print("---Output GMT---")


def main():
    BACKUP_DIR.mkdir(exist_ok=True)

    paths_to_restore = [
        DATA_DIR / "input.py",
        DATA_DIR / "Input_soil_material_parameters.py",
        SRC_DIR / "input",
        RESULTS_DIR,
    ]

    try:
        for path in paths_to_restore:
            backup_path(path)

        copy_case_files()
        run_model()
        run_gmt_plots()
        collect_results()
        plotting()

        print(f"Output saved to: {BENCH_OUTPUT_DIR}")

    finally:
        for path in paths_to_restore:
            restore_path(path)

        shutil.rmtree(BACKUP_DIR, ignore_errors=True)


if __name__ == "__main__":
    main()