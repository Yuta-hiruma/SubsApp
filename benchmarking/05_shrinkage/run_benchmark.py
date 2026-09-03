from pathlib import Path
import shutil
import subprocess
import sys


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
    Se_vG_fig = RESULTS_DIR /  "figure/Se_vG.png"
    #benchmarking---------------------
    stress_strain_fig = RESULTS_DIR /  "figure/stress_strain.png"
    pressure_void_fig = RESULTS_DIR /  "figure/pressure_void.png"
    pressure_strain_fig = RESULTS_DIR /  "figure/pressure_strain.png"
    void_time_fig = RESULTS_DIR /  "figure/void_time.png"
    pressure_time_fig = RESULTS_DIR /  "figure/pressure_time.png"
    sscc_fig = RESULTS_DIR /  "figure/sscc.png"
    sscc_vG_fig = RESULTS_DIR /  "figure/sscc_vG.png"
    sscc_MR_fig = RESULTS_DIR /  "figure/sscc_MR.png"
    #gwl_fig = RESULTS_DIR /  "figure/gwl.png"


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
    shutil.copy2(stress_strain_fig, BENCH_OUTPUT_DIR / "stress_strain.png",)
    shutil.copy2(pressure_strain_fig, BENCH_OUTPUT_DIR / "pressure_strain.png",)
    shutil.copy2(pressure_void_fig, BENCH_OUTPUT_DIR / "pressure_void.png",)
    shutil.copy2(void_time_fig, BENCH_OUTPUT_DIR / "void_time.png",)
    shutil.copy2(pressure_time_fig, BENCH_OUTPUT_DIR / "pressure_time.png",)
    shutil.copy2(sscc_fig, BENCH_OUTPUT_DIR / "sscc.png",)
    shutil.copy2(sscc_vG_fig, BENCH_OUTPUT_DIR / "sscc_vG.png",)
    shutil.copy2(sscc_MR_fig, BENCH_OUTPUT_DIR / "sscc_MR.png",)
    #shutil.copy2(gwl_fig, BENCH_OUTPUT_DIR / "gwl.png",)


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

        print(f"Output saved to: {BENCH_OUTPUT_DIR}")

    finally:
        for path in paths_to_restore:
            restore_path(path)

        shutil.rmtree(BACKUP_DIR, ignore_errors=True)


if __name__ == "__main__":
    main()