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
    shutil.copy2(displacement_file, BENCH_OUTPUT_DIR / "displacement.dat",)



def plotting():
    # read measurement data-----
    measure_displacement = "sub-cr_result/sub-cr_ch5.csv"
    displacement_bp = pd.read_csv(measure_displacement, sep=",", skiprows=3, names=["time","ocr15","ocr17","ocr20"])
    
    # read simulation result-----
    displacement_file = RESULTS_DIR / "displacement.dat"
    displacement = pd.read_csv(displacement_file,sep=r"\s+",skiprows=1,names=["days", "displacement"])
    displacement["time"] = pd.to_datetime(displacement["days"])
    displacement["time"] = (displacement["time"] - displacement["time"].iloc[0]).dt.days

    #-----------------------------------
    #displacement
    plt.figure()
    plt.plot(displacement_bp["time"],displacement_bp["ocr20"],marker="s",linestyle=":",label="Sub-CR")
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
    plt.savefig(BENCH_OUTPUT_DIR / "displacement.png")#,dpi=300)
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