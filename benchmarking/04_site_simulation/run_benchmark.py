from pathlib import Path
import shutil
import subprocess
import sys
import pandas as pd
import matplotlib.pyplot as plt

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

    Se_vG_fig = RESULTS_DIR /  "figure/Se_vG.png"
    stress_strain_fig = RESULTS_DIR /  "figure/stress_strain.png"
    sscc_MR_fig = RESULTS_DIR /  "figure/sscc_MR.png"
    Pressure_fig = RESULTS_DIR /  "figure/pressure.png"
    theta_file = RESULTS_DIR / "theta.dat"
    theta_fig = RESULTS_DIR /  "figure/theta.png"



    shutil.copy2(Se_vG_fig, BENCH_OUTPUT_DIR / "Se_vG.png",)
    shutil.copy2(stress_strain_fig, BENCH_OUTPUT_DIR / "stress_strain.png",)
    shutil.copy2(sscc_MR_fig, BENCH_OUTPUT_DIR / "sscc_MR.png",)
    shutil.copy2(Pressure_fig, BENCH_OUTPUT_DIR / "pressure.png",)
    shutil.copy2(theta_file, BENCH_OUTPUT_DIR / "theta.dat",)
    shutil.copy2(theta_fig, BENCH_OUTPUT_DIR / "theta.png",)

def plotting():
    #read measurement -----------
    measured_displacement_file = "measurement/lvl.dat"
    measured_displacement = pd.read_csv(measured_displacement_file,sep=r"\s+",skiprows=1,names=["date", "EXT-K1", "EXT-K2", "EXT-K3", "EXT-K4", "EXT-K5", "EXT-K6"])
    measured_displacement["EXT-K6"] = (measured_displacement["EXT-K6"] - measured_displacement["EXT-K6"].iloc[0])
    # read simulation result-----
    stress_strain_file = RESULTS_DIR / "stress_strain.dat"
    displacement_file = RESULTS_DIR / "displacement.dat"
    sscc_MR_file = RESULTS_DIR / "sscc_MR_simulation.dat"
    stress_strain = pd.read_csv(stress_strain_file,sep=r"\s+",skiprows=1,names=["sigma_eff", "strain"])
    displacement = pd.read_csv(displacement_file,sep=r"\s+",skiprows=1,names=["time", "displacement"])
    sscc_MR = pd.read_csv(sscc_MR_file,sep=r"\s+",skiprows=1,names=["moisture_ratio", "void"])

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
    plt.plot(displacement["time"].iloc[:-1],measured_displacement["EXT-K6"]*10,label="Measurement")
    plt.title("Vertical displacement at surface ")
    plt.xlabel(r"time [day]")
    plt.ylabel(r"Vertical displacement [$cm$]")
    #plt.grid(True,which="both",linestyle="--",alpha=0.5)
    plt.legend()
    plt.tight_layout()
    plt.savefig(BENCH_OUTPUT_DIR / "displacement.png",dpi=300)
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
    plt.savefig(BENCH_OUTPUT_DIR / "sscc_MR.png",dpi=300)
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