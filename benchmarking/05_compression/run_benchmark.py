from pathlib import Path
import shutil
import subprocess
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

#Read input------------
current_dir = os.getcwd()
parent_dir = os.path.abspath("data")
sys.path.insert(0, parent_dir)
from Input_soil_material_parameters import *
from input import *
sys.path.pop(0)
os.chdir(current_dir)
pressure_ini_file =  "input/pressure_ini/pressure_ini_swp.dat"
pressure_ini = pd.read_csv(pressure_ini_file,sep=r"\s+",skiprows=1,names=["depth", "pressure"])
sigma_eff_0 = pressure_ini["pressure"][0]
OCR = Soil_parameters["Clay"]["Compression"]["OCR"]
sigma_eff_p = sigma_eff_0 * OCR
Kappa =Soil_parameters["Clay"]["Compression"]["kappa"]
Lambda = Soil_parameters["Clay"]["Compression"]["lambda"]
Mu = Soil_parameters["Clay"]["Compression"]["mu"]
print("sigma_eff_0 =",sigma_eff_0, "kPa")
print("sigma_eff_p =",sigma_eff_p, "kPa")

#------------------------

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
    Pressure_fig = RESULTS_DIR /  "gmt/pressure.png"
    theta_file = RESULTS_DIR / "theta.dat"
    theta_fig = RESULTS_DIR /  "gmt/theta.png"
    Se_vG_fig = RESULTS_DIR /  "gmt/Se_vG.png"
    #benchmarking---------------------
    #stress_strain_fig = RESULTS_DIR /  "figure/stress_strain.png"
    pressure_void_fig = RESULTS_DIR /  "figure/pressure_void.png"
    pressure_time_fig = RESULTS_DIR /  "figure/pressure_time.png"
    sscc_MR_fig = RESULTS_DIR /  "figure/sscc_MR.png"

    if not displacement_file.exists():
        raise FileNotFoundError(
            f"Expected output was not created: {displacement_file}"
        )

    shutil.copy2(displacement_file, BENCH_OUTPUT_DIR / "displacement.dat",)
    shutil.copy2(displacement_fig, BENCH_OUTPUT_DIR / "displacement.png",)
    shutil.copy2(Pressure_top_file, BENCH_OUTPUT_DIR / "Pressure_top.dat",)
    shutil.copy2(Pressure_bottom_file, BENCH_OUTPUT_DIR / "Pressure_bottom.dat",)
    shutil.copy2(Pressure_fig, BENCH_OUTPUT_DIR / "pressure.png",)
    shutil.copy2(theta_file, BENCH_OUTPUT_DIR / "theta.dat",)
    shutil.copy2(theta_fig, BENCH_OUTPUT_DIR / "theta.png",)
    shutil.copy2(Se_vG_fig, BENCH_OUTPUT_DIR / "Se_vG.png",)
    #benchmarking---------------------
    #shutil.copy2(stress_strain_fig, BENCH_OUTPUT_DIR / "stress_strain.png",)
    shutil.copy2(pressure_void_fig, BENCH_OUTPUT_DIR / "pressure_void.png",)
    shutil.copy2(pressure_time_fig, BENCH_OUTPUT_DIR / "pressure_time.png",)
    shutil.copy2(sscc_MR_fig, BENCH_OUTPUT_DIR / "sscc_MR.png",)


def plotting():
    # read simulation result-----
    stress_strain_file = RESULTS_DIR / "stress_strain.dat"
    displacement_file = RESULTS_DIR / "displacement.dat"
    sscc_MR_file = RESULTS_DIR / "sscc_MR_simulation.dat"
    
    stress_strain = pd.read_csv(stress_strain_file,sep=r"\s+",skiprows=1,names=["sigma_eff", "strain", "void_ratio"])
    displacement = pd.read_csv(displacement_file,sep=r"\s+",skiprows=1,names=["time", "displacement"])
    sscc_MR = pd.read_csv(sscc_MR_file,sep=r"\s+",skiprows=1,names=["moisture_ratio", "void"])

    #-----------------------------------
    #Effective stress - strain
    sigma = stress_strain["sigma_eff"].to_numpy()
    strain = -stress_strain["strain"].to_numpy()
    void_ratio = stress_strain["void_ratio"].to_numpy()
    order = np.argsort(sigma)
    sigma_sorted = sigma[order]
    strain_sorted = strain[order]

    
    strain_0 = np.interp(0, sigma_sorted, strain_sorted)
    
    # Swelling/recompression line
    stress_kappa = np.geomspace(sigma_sorted.min(), sigma_eff_p, 100)
    void_kappa = np.full(len(stress_kappa), void_ratio[0])
    strain_kappa = Kappa/(1+void_kappa) * np.log(stress_kappa / sigma_eff_0)
    strain_p = np.max(strain_kappa)
    
    # Virgin compression line
    stress_lambda = np.geomspace(sigma_eff_p, sigma_sorted.max(), 100)
    void_lambda = np.full(len(stress_lambda), void_ratio[-1])
    strain_lambda = strain_p + Lambda/(1+void_lambda) * np.log(stress_lambda / sigma_eff_p)


    fig, ax = plt.subplots()
    ax.plot(stress_strain["sigma_eff"],-stress_strain["strain"],ls="-",color="black",label="Simulation",linewidth=0.8,)
    ax.invert_yaxis()
    ax.axvline(x=sigma_eff_0,color="black",linestyle=":")
    ax.axvline(x=sigma_eff_p,color="black",linestyle=":")
    ax.text(sigma_eff_0, 0.03,r"$\sigma'_{0}$",transform=ax.get_xaxis_transform(),ha="left", va="bottom", color="black",)
    ax.text(sigma_eff_p, 0.03,r"$\sigma'_{p}$",transform=ax.get_xaxis_transform(),ha="left", va="bottom", color="black",)
    ax.plot(stress_kappa, strain_kappa, "--", color="blue", alpha=0.6)
    ax.plot(stress_lambda, strain_lambda, "--", color="orange", alpha=0.6)
    i = len(stress_kappa) // 2
    ax.annotate(r"$\kappa^*$",(stress_kappa[i], strain_kappa[i]),xytext=(0, 8), textcoords="offset points",color="blue", ha="center",)
    i = len(stress_lambda) // 2
    ax.annotate(r"$\lambda^*$",(stress_lambda[i], strain_lambda[i]),xytext=(0, 8), textcoords="offset points",color="orange", ha="center",)
    #ax.set_title("Stress–strain curve")
    ax.set_xlabel(r"Effective stress $\sigma'$ [kPa]")
    ax.set_ylabel(r"Strain $\varepsilon$ [-]")
    ax.set_xscale("log")
    ax.grid(True, which="both", linestyle="--", alpha=0.5)
    fig.tight_layout()
    fig.savefig(BENCH_OUTPUT_DIR / "stress_strain.png")
    plt.close(fig)
    
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
    plt.plot(sscc_MR["moisture_ratio"],sscc_MR["void"],label="Simulation")
    plt.plot(sscc_MR["void"],sscc_MR["void"],linestyle="--",color="black",label="Saturated soil")
    plt.title("Soil shrinkage characteristic curve")
    plt.xlabel(r"Moisture ratio $\vartheta$")
    plt.ylabel(r"Void ratio $e$")
    plt.legend()
    plt.grid(True,which="both",linestyle="--",alpha=0.5)
    plt.tight_layout()
    plt.savefig(BENCH_OUTPUT_DIR / "sscc_MR.png")
    plt.close()


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
        collect_results()
        plotting()

        print(f"Output saved to: {BENCH_OUTPUT_DIR}")

    finally:
        for path in paths_to_restore:
            restore_path(path)

        shutil.rmtree(BACKUP_DIR, ignore_errors=True)


if __name__ == "__main__":
    main()