from pathlib import Path
import shutil
import subprocess
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
current_dir = os.getcwd()
parent_dir = os.path.abspath("data")
sys.path.insert(0, parent_dir)
from Input_soil_material_parameters import *
from input import *
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
    Se_vG_fig = RESULTS_DIR /  "figure/Se_vG.png"
    stress_strain_fig = RESULTS_DIR /  "figure/stress_strain.png"
    sscc_MR_fig = RESULTS_DIR /  "figure/sscc_MR.png"
    Pressure_fig = RESULTS_DIR /  "figure/pressure.png"
    Pressure_file = RESULTS_DIR /  "Pressure_top.dat"
    theta_file = RESULTS_DIR / "theta.dat"
    theta_fig = RESULTS_DIR /  "figure/theta.png"

    shutil.copy2(displacement_file, BENCH_OUTPUT_DIR / "displacement.dat",)
    shutil.copy2(Se_vG_fig, BENCH_OUTPUT_DIR / "Se_vG.png",)
    shutil.copy2(stress_strain_fig, BENCH_OUTPUT_DIR / "stress_strain.png",)
    shutil.copy2(sscc_MR_fig, BENCH_OUTPUT_DIR / "sscc_MR.png",)
    shutil.copy2(Pressure_fig, BENCH_OUTPUT_DIR / "pressure.png",)
    shutil.copy2(Pressure_file, BENCH_OUTPUT_DIR / "Pressure_top.dat",)
    shutil.copy2(theta_file, BENCH_OUTPUT_DIR / "theta.dat",)
    shutil.copy2(theta_fig, BENCH_OUTPUT_DIR / "theta.png",)

def plotting():
    #read measurement -----------
    measured_displacement_file = "measurement/lvl.dat"
    measured_gwl_file = "measurement/groundwater_level.dat"
    measured_displacement = pd.read_csv(measured_displacement_file,sep=r"\s+",skiprows=1,names=["date", "EXT-K1", "EXT-K2", "EXT-K3", "EXT-K4", "EXT-K5", "EXT-K6"])
    measured_displacement["EXT-K1_0"] = (measured_displacement["EXT-K1"] - measured_displacement["EXT-K1"].iloc[0])
    measured_displacement["EXT-K2_0"] = (measured_displacement["EXT-K2"] - measured_displacement["EXT-K2"].iloc[0])
    measured_displacement["EXT-K3_0"] = (measured_displacement["EXT-K3"] - measured_displacement["EXT-K3"].iloc[0])
    measured_displacement["EXT-K4_0"] = (measured_displacement["EXT-K4"] - measured_displacement["EXT-K4"].iloc[0])
    measured_displacement["EXT-K5_0"] = (measured_displacement["EXT-K5"] - measured_displacement["EXT-K5"].iloc[0])
    measured_displacement["surface"] = (measured_displacement["EXT-K6"] - measured_displacement["EXT-K6"].iloc[0])
    measured_gwl = pd.read_csv(measured_gwl_file,sep=r"\s+",skiprows=1,names=["date", "gwl"])
    measured_gwl["time"] = pd.to_datetime(measured_gwl["date"])
    measured_gwl["time"] = (measured_gwl["time"] - measured_gwl["time"].iloc[0]).dt.days

    # read simulation result-----
    stress_strain_file = RESULTS_DIR / "stress_strain.dat"
    displacement_file = RESULTS_DIR / "displacement.dat"
    node_heights_file = RESULTS_DIR / "node_heights.dat"
    sscc_MR_file = RESULTS_DIR / "sscc_MR_simulation.dat"
    gwl_file = RESULTS_DIR / "gwl.dat"
    stress_strain = pd.read_csv(stress_strain_file,sep=r"\s+",skiprows=1,names=["sigma_eff", "strain"])
    displacement = pd.read_csv(displacement_file,sep=r"\s+",skiprows=1,names=["time", "displacement"])
    sscc_MR = pd.read_csv(sscc_MR_file,sep=r"\s+",skiprows=1,names=["moisture_ratio", "void"])
    gwl = pd.read_csv(gwl_file,sep=r"\s+",skiprows=1,names=["time", "gwl"])

    node_heights = pd.read_csv(node_heights_file,sep=r"\s+",skiprows=1)
    node_time = node_heights.iloc[:, 0].to_numpy()
    node_height = node_heights.iloc[:, 1:].to_numpy()



    #-----------------------------------
    #Effective stress - strain
    plt.figure()
    plt.plot(stress_strain["sigma_eff"],stress_strain["strain"],".",color="black")
    plt.title("Stress strain")
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
    plt.plot(displacement["time"].iloc[:-1],measured_displacement["surface"],color="black",label="Measurement")
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
    plt.title("Soil shrinkage characteristic curve")
    plt.xlabel(r"Moisture ratio $\vartheta$")
    plt.ylabel(r"Void ratio $e$")
    plt.legend()
    plt.grid(True,which="both",linestyle="--",alpha=0.5)

    plt.tight_layout()
    plt.savefig(BENCH_OUTPUT_DIR / "sscc_MR.png")
    plt.close()

    #-----------------------------------
    #Groundwater level
    plt.figure()
    plt.plot(gwl["time"],gwl["gwl"],label="Simulation")
    plt.plot(measured_gwl["time"],measured_gwl["gwl"]/100,label="Measurement",linestyle="--",color="blue")
    plt.title("Groundwater level")
    plt.xlabel(r"Time [day]")
    plt.ylabel(r"Groundwater level[m]")
    plt.legend()
    plt.grid(True,which="both",linestyle="--",alpha=0.5)
    plt.tight_layout()
    plt.savefig(BENCH_OUTPUT_DIR / "gwl.png")
    plt.close()

    ext1 = abs(int(measured_displacement["EXT-K1"].iloc[0]/100))
    ext2 = abs(int(measured_displacement["EXT-K2"].iloc[0]/100))
    ext3 = abs(int(measured_displacement["EXT-K3"].iloc[0]/100))
    ext4 = abs(int(measured_displacement["EXT-K4"].iloc[0]/100))
    ext5 = abs(int(measured_displacement["EXT-K5"].iloc[0]/100))
 

    #-----------------------------------
    #Vertical displacement at extensometer
    plt.figure()
    plt.plot(node_time,np.sum(node_height[:, ext5:-1]-node_height[0, ext5:-1], axis=1) ,label="Simulation")
    plt.plot(displacement["time"].iloc[:-1],measured_displacement["EXT-K5_0"],color="black",label="Measurement")
    plt.title("Vertical displacement at surface ")
    plt.xlabel(r"time [day]")
    plt.ylabel(r"Vertical displacement at EXT-K5 [$mm$]")
    #plt.grid(True,which="both",linestyle="--",alpha=0.5)
    plt.legend()
    plt.subplots_adjust(left=0.2, bottom=0.2)
    plt.gca().xaxis.set_major_locator(plt.MaxNLocator(8))
    plt.tight_layout()
    plt.savefig(BENCH_OUTPUT_DIR / "displacement_ext5.png")#,dpi=300)
    plt.close()  


    plt.figure()
    plt.plot(node_time,np.sum(node_height[:, ext4:-1]-node_height[0, ext4:-1], axis=1) ,label="Simulation")
    plt.plot(displacement["time"].iloc[:-1],measured_displacement["EXT-K4_0"],color="black",label="Measurement")
    plt.title("Vertical displacement at surface ")
    plt.xlabel(r"time [day]")
    plt.ylabel(r"Vertical displacement at EXT-K4 [$mm$]")
    #plt.grid(True,which="both",linestyle="--",alpha=0.5)
    plt.legend()
    plt.subplots_adjust(left=0.2, bottom=0.2)
    plt.gca().xaxis.set_major_locator(plt.MaxNLocator(8))
    plt.tight_layout()
    plt.savefig(BENCH_OUTPUT_DIR / "displacement_ext4.png")#,dpi=300)
    plt.close()  

    plt.figure()
    plt.plot(node_time,np.sum(node_height[:, ext3:-1]-node_height[0, ext3:-1], axis=1) ,label="Simulation")
    plt.plot(displacement["time"].iloc[:-1],measured_displacement["EXT-K3_0"],color="black",label="Measurement")
    plt.title("Vertical displacement at surface ")
    plt.xlabel(r"time [day]")
    plt.ylabel(r"Vertical displacement at EXT-K3 [$mm$]")
    #plt.grid(True,which="both",linestyle="--",alpha=0.5)
    plt.legend()
    plt.subplots_adjust(left=0.2, bottom=0.2)
    plt.gca().xaxis.set_major_locator(plt.MaxNLocator(8))
    plt.tight_layout()
    plt.savefig(BENCH_OUTPUT_DIR / "displacement_ext3.png")#,dpi=300)
    plt.close()  

    plt.figure()
    plt.plot(node_time,np.sum(node_height[:, ext2:-1]-node_height[0, ext2:-1], axis=1) ,label="Simulation")
    plt.plot(displacement["time"].iloc[:-1],measured_displacement["EXT-K2_0"],color="black",label="Measurement")
    plt.title("Vertical displacement at surface ")
    plt.xlabel(r"time [day]")
    plt.ylabel(r"Vertical displacement at EXT-K2 [$mm$]")
    #plt.grid(True,which="both",linestyle="--",alpha=0.5)
    plt.legend()
    plt.subplots_adjust(left=0.2, bottom=0.2)
    plt.gca().xaxis.set_major_locator(plt.MaxNLocator(8))
    plt.tight_layout()
    plt.savefig(BENCH_OUTPUT_DIR / "displacement_ext2.png")#,dpi=300)
    plt.close()  

    plt.figure()
    plt.plot(node_time,np.sum(node_height[:, ext1:-1]-node_height[0, ext1:-1], axis=1) ,label="Simulation")
    plt.plot(displacement["time"].iloc[:-1],measured_displacement["EXT-K1_0"],color="black",label="Measurement")
    plt.title("Vertical displacement at surface ")
    plt.xlabel(r"time [day]")
    plt.ylabel(r"Vertical displacement at EXT-K1 [$mm$]")
    #plt.grid(True,which="both",linestyle="--",alpha=0.5)
    plt.legend()
    plt.subplots_adjust(left=0.2, bottom=0.2)
    plt.gca().xaxis.set_major_locator(plt.MaxNLocator(8))
    plt.tight_layout()
    plt.savefig(BENCH_OUTPUT_DIR / "displacement_ext1.png")#,dpi=300)
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