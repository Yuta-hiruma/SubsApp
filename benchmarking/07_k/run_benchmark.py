from pathlib import Path
import shutil
import subprocess
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib import ticker

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

    displacement_file = RESULTS_DIR / "displacement.dat"
    displacement_fig = RESULTS_DIR /  "figure/displacement.png"
    Pressure_top_file = RESULTS_DIR / "Pressure_top.dat"
    Pressure_bottom_file = RESULTS_DIR / "Pressure_bottom.dat"
    theta_file = RESULTS_DIR / "theta.dat"
    k_file = RESULTS_DIR / "k.dat"
    pressure_fig = RESULTS_DIR /  "gmt/pressure.png"
    theta_fig = RESULTS_DIR /  "gmt/theta.png"
    Se_vG_fig = RESULTS_DIR /  "gmt/Se_vG.png"
    K_fig = RESULTS_DIR /  "gmt/K.png"
    Se_vG_file = RESULTS_DIR / "Se_vG.dat"
    K_file = RESULTS_DIR / "K.dat"


    if not displacement_file.exists():
        raise FileNotFoundError(
            f"Expected output was not created: {displacement_file}"
        )

    shutil.copy2(displacement_file, BENCH_OUTPUT_DIR / "displacement.dat",)
    shutil.copy2(displacement_fig, BENCH_OUTPUT_DIR / "displacement.png",)
    shutil.copy2(Pressure_top_file, BENCH_OUTPUT_DIR / "Pressure_top.dat",)
    shutil.copy2(Pressure_bottom_file, BENCH_OUTPUT_DIR / "Pressure_bottom.dat",)
    shutil.copy2(pressure_fig, BENCH_OUTPUT_DIR / "pressure.png",)
    shutil.copy2(theta_file, BENCH_OUTPUT_DIR / "theta.dat",)
    shutil.copy2(k_file, BENCH_OUTPUT_DIR / "k.dat",)
    shutil.copy2(theta_fig, BENCH_OUTPUT_DIR / "theta.png",)
    shutil.copy2(Se_vG_fig, BENCH_OUTPUT_DIR / "Se_vG.png",)
    shutil.copy2(Se_vG_file, BENCH_OUTPUT_DIR / "Se_vG.dat",)
    shutil.copy2(K_fig, BENCH_OUTPUT_DIR / "K.png",)
    shutil.copy2(K_file, BENCH_OUTPUT_DIR / "K.dat",)

def plotting():
    # read simulation result-----
    stress_strain_file = RESULTS_DIR / "stress_strain.dat"
    displacement_file = RESULTS_DIR / "displacement.dat"
    sscc_MR_file = RESULTS_DIR / "sscc_MR_simulation.dat"
    Se_vG_file = RESULTS_DIR / "Se_vG.dat"
    stress_strain = pd.read_csv(stress_strain_file,sep=r"\s+",skiprows=1,names=["sigma_eff", "strain"])
    displacement = pd.read_csv(displacement_file,sep=r"\s+",skiprows=1,names=["time", "displacement"])
    sscc_MR = pd.read_csv(sscc_MR_file,sep=r"\s+",skiprows=1,names=["moisture_ratio", "void"])
    Se_vG = pd.read_csv(Se_vG_file,sep=r"\s+",skiprows=1,names=["time", "depth", "effective_saturation"],)
    Se_grid = Se_vG.pivot(index="depth",columns="time",values="effective_saturation",).sort_index().sort_index(axis=1)

    # -----------------------------------
    # Se_vG
    Se_vG_file = RESULTS_DIR / "Se_vG.dat"
    Se_vG = pd.read_csv(Se_vG_file,sep=r"\s+",skiprows=1,names=["time", "depth", "effective_saturation"],)
    Se_grid = Se_vG.pivot(index="depth",columns="time",values="effective_saturation",).sort_index().sort_index(axis=1)
    simulation_time = Se_grid.columns.to_numpy(dtype=float)
    depth = Se_grid.index.to_numpy(dtype=float)
    Se_values = Se_grid.to_numpy(dtype=float)
 
    Se_min = Se_values.min()
    Se_max = Se_values.max()
    
    bounds = np.linspace(Se_min, Se_max, 31)
    if Se_min >= 1.0:
        bounds = np.linspace(0.0, max(1.0, Se_max), 31)
    
    colors = plt.get_cmap("Blues")(np.linspace(0, 1, 30))
    colors[0] = (1, 1, 1, 1)
    cmap = ListedColormap(colors)
    norm = BoundaryNorm(bounds, ncolors=cmap.N, clip=True)

    plt.pcolormesh(simulation_time,depth,Se_values, shading="nearest",cmap=cmap, norm=norm)

    #-----------------------------------
    #Effective stress - strain
    plt.figure()
    plt.plot(stress_strain["sigma_eff"],stress_strain["strain"],".",color="black")
    #plt.title("Stress strain")
    plt.xlabel(r"Effective stress $\sigma'$ [kPa]")
    plt.ylabel(r"Strain $\varepsilon$ [-]")
    plt.xscale("log")
    plt.grid(True,which="both",linestyle="--",alpha=0.5)
    plt.tight_layout()
    plt.savefig(BENCH_OUTPUT_DIR / "stress_strain.png",dpi=300)
    plt.close()

    #-----------------------------------
    #Vertical displacement at surface
    plt.figure()
    plt.plot(displacement["time"],displacement["displacement"],label="Simulation")
    plt.title("Vertical displacement at surface ")
    plt.xlabel(r"time [day]")
    plt.ylabel(r"Vertical displacement [$mm$]")
    #plt.grid(True,which="both",linestyle="--",alpha=0.5)
    plt.legend()
    plt.subplots_adjust(left=0.2, bottom=0.2)
    plt.gca().xaxis.set_major_locator(plt.MaxNLocator(8))
    plt.tight_layout()
    plt.savefig(BENCH_OUTPUT_DIR / "displacement.png")#,dpi=300)
    plt.close()  

    #-----------------------------------
    #Soil shrinkage characteristic curve
    plt.figure()
    plt.plot(sscc_MR["moisture_ratio"],sscc_MR["void"],"." ,color = "black",label="Simulation" )
    plt.plot(sscc_MR["void"],sscc_MR["void"],linestyle="--",color="black",label="Saturated soil")
    #plt.title("Soil shrinkage characteristic curve")
    plt.xlabel(r"Moisture ratio $\vartheta [-]$")
    plt.ylabel(r"Void ratio $e[-]$")
    plt.legend()
    plt.grid(True,which="both",linestyle="--",alpha=0.5)
    plt.tight_layout()
    plt.savefig(BENCH_OUTPUT_DIR / "sscc_MR.png")
    plt.close()

def run_gmt_plots():
    commands = [
        ["python3", "plotgmt.py"],
        ["bash", "plot_Se_vG.sh"],
        ["bash", "plot_K.sh"],
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