from pathlib import Path
import shutil
import subprocess
import sys


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
    pressure_fig = RESULTS_DIR /  "gmt/pressure.png"
    theta_fig = RESULTS_DIR /  "gmt/theta.png"
    Se_vG_fig = RESULTS_DIR /  "gmt/Se_vG.png"


    if not displacement_file.exists():
        raise FileNotFoundError(
            f"Expected output was not created: {displacement_file}"
        )

    shutil.copy2(displacement_file, BENCH_OUTPUT_DIR / "displacement.dat",)
    shutil.copy2(displacement_fig, BENCH_OUTPUT_DIR / "displacement.png",)
    shutil.copy2(Pressure_top_file, BENCH_OUTPUT_DIR / "Pressure_top.dat",)
    shutil.copy2(Pressure_bottom_file, BENCH_OUTPUT_DIR / "Pressure_bottom.dat",)
    shutil.copy2(theta_file, BENCH_OUTPUT_DIR / "theta.dat",)
    shutil.copy2(pressure_fig, BENCH_OUTPUT_DIR / "pressure.png",)
    shutil.copy2(theta_fig, BENCH_OUTPUT_DIR / "theta.png",)
    shutil.copy2(Se_vG_fig, BENCH_OUTPUT_DIR / "Se_vG.png",)

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

        print(f"Output saved to: {BENCH_OUTPUT_DIR}")

    finally:
        for path in paths_to_restore:
            restore_path(path)

        shutil.rmtree(BACKUP_DIR, ignore_errors=True)


if __name__ == "__main__":
    main()