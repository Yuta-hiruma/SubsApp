from pathlib import Path
import shutil
import subprocess
import sys
import pandas as pd
import matplotlib.pyplot as plt
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
    Pressure_fig = RESULTS_DIR /  "figure/pressure.png"
    theta_file = RESULTS_DIR / "theta.dat"
    theta_fig = RESULTS_DIR /  "figure/theta.png"

    shutil.copy2(displacement_file, BENCH_OUTPUT_DIR / "displacement.dat",)
    shutil.copy2(displacement_fig, BENCH_OUTPUT_DIR / "displacement.png",)
    shutil.copy2(Pressure_top_file, BENCH_OUTPUT_DIR / "Pressure_top.dat",)
    shutil.copy2(Pressure_bottom_file, BENCH_OUTPUT_DIR / "Pressure_bottom.dat",)
    shutil.copy2(Pressure_fig, BENCH_OUTPUT_DIR / "pressure.png",)
    shutil.copy2(theta_file, BENCH_OUTPUT_DIR / "theta.dat",)
    shutil.copy2(theta_fig, BENCH_OUTPUT_DIR / "theta.png",)

def plotting():
    # read measurement data-----
    measure_displacement = "sub-cr_result/sub-cr_ch5.csv"
    displacement_bp = pd.read_csv(measure_displacement, sep=",", skiprows=3, names=["time","ocr15","ocr17","ocr20"])
    
    # read simulation result-----
    displacement_file = RESULTS_DIR / "displacement.dat"
    displacement = pd.read_csv(displacement_file,sep=r"\s+",skiprows=1,names=["days", "displacement"])
    displacement["time"] = pd.to_datetime(displacement["days"])
    displacement["time"] = (displacement["time"] - displacement["time"].iloc[0]).dt.days
    stress_strain_file = RESULTS_DIR / "stress_strain.dat"
    pressure_strain_file = RESULTS_DIR / "pressure_strain_simulation.dat"
    pressure_void_file = RESULTS_DIR / "pressure_void_simulation.dat"
    void_time_file = RESULTS_DIR / "void_time_simulation.dat"
    pressure_time_file = RESULTS_DIR / "pressure_time_simulation.dat"
    stress_strain = pd.read_csv(stress_strain_file,sep=r"\s+",skiprows=1,names=["sigma_eff", "strain"])
    pressure_strain = pd.read_csv(pressure_strain_file,sep=r"\s+",skiprows=1,names=["pressure", "strain"])
    pressure_void = pd.read_csv(pressure_void_file,sep=r"\s+",skiprows=1,names=["pressure", "void"])
    void_time = pd.read_csv(void_time_file,sep=r"\s+",skiprows=1,names=["time", "void"])
    pressure_time = pd.read_csv(pressure_time_file,sep=r"\s+",skiprows=1,names=["time", "pressure"])

    #-----------------------------------
    #displacement
    plt.figure()
    plt.plot(displacement_bp["time"],displacement_bp["ocr15"],marker="s",linestyle=":",label="Sub-CR")
    plt.plot(displacement["time"],displacement["displacement"]/1000,label="SubsApp")
    plt.title("Comparison between Sub-CR and SubsApp")
    plt.xlabel(r"Time $\psi$ [day]")
    plt.ylabel(r"Vetical displacement [m]")
    plt.legend()
    plt.grid(True,which="both",linestyle="--",alpha=0.5)
    plt.tight_layout()
    textstr = (
    rf"$OCR = {Soil_parameters["confining_layer"]["Compression"]["OCR"]:.3f}$" "\n"
    )
    plt.text(0.05, 0.05,textstr,transform=plt.gca().transAxes,fontsize=11,verticalalignment="bottom",
    bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))
    plt.savefig(BENCH_OUTPUT_DIR / "displacement.png",dpi=300)
    plt.close()

    #-----------------------------------
    #Effective stress - strain
    plt.figure()
    plt.plot(stress_strain["sigma_eff"],stress_strain["strain"])#,marker="s",linestyle=":")
    plt.title("Stress strain")
    plt.xlabel(r"Effective stress $\sigma'$ [kPa]")
    plt.ylabel(r"Strain $\varepsilon$ [-]")
    plt.xticks(rotation=45)
    plt.xscale("log")
    plt.grid(True,which="both",linestyle="--",alpha=1.0)
    plt.tight_layout()

    textstr = (
    rf"$\kappa = {Soil_parameters["confining_layer"]["Compression"]["kappa"]:.3f}$" "\n"
    rf"$\lambda = {Soil_parameters["confining_layer"]["Compression"]["lambda"]:.3f}$" "\n"
    rf"$\mu = {Soil_parameters["confining_layer"]["Compression"]["mu"]:.3f}$"
    )
    #plt.text(0.05, 0.05,textstr,transform=plt.gca().transAxes,fontsize=11,verticalalignment="bottom",bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))
    plt.savefig(BENCH_OUTPUT_DIR / "stress_strain.png",dpi=300)
    plt.close()

    #-----------------------------------
    #Pore pressure - void ratio
    plt.figure()
    plt.plot(pressure_void["pressure"],pressure_void["void"],marker="s",linestyle=":")
    plt.title("Pressure - void ratio")
    plt.xlabel(r"Pore pressure $\psi$ [kPa]")
    plt.ylabel(r"Void ratio $e$ [-]")
    #plt.legend()
    plt.grid(True,which="both",linestyle="--",alpha=0.5)
    plt.tight_layout()
    plt.savefig(BENCH_OUTPUT_DIR / "pressure_void.png",dpi=300)
    plt.close()

    #-----------------------------------
    #Pore pressure - strain
    plt.figure()
    plt.plot(pressure_strain["pressure"],pressure_strain["strain"],marker="s",linestyle=":")
    plt.title("Pressure - strain")
    plt.xlabel(r"Pore pressure $\psi$ [kPa]")
    plt.ylabel(r"Strain $\varepsilon$ [-]")
    plt.grid(True,which="both",linestyle="--",alpha=0.5)
    plt.tight_layout()
    plt.savefig(BENCH_OUTPUT_DIR / "pressure_strain.png",dpi=300)
    plt.close()

    #-----------------------------------
    #Void ratio - time
    plt.figure()
    plt.plot(void_time["time"],void_time["void"])
    plt.title("Void ratio - time")
    plt.xlabel("Time [day]")
    plt.ylabel(r"Void ratio $e$")
    plt.grid(True,which="both",linestyle="--",alpha=0.5)
    plt.tight_layout()
    plt.savefig(BENCH_OUTPUT_DIR / "void_time.png",dpi=300)
    plt.close()
    #-----------------------------------
    #Pore pressure - time
    plt.figure()
    plt.plot(pressure_time["time"],pressure_time["pressure"],marker="s",linestyle=":")
    plt.title("Pressure - time")
    plt.xlabel("Time [day]")
    plt.ylabel("Pore pressure [kPa]")
    plt.grid(True,which="both",linestyle="--",alpha=0.5)
    plt.tight_layout()
    plt.savefig(BENCH_OUTPUT_DIR / "pressure_time.png",dpi=300)
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