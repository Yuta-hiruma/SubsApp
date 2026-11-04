import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr

#read measurement -----------
measured_displacement_file = "measurement/lvl.dat"
measured_displacement = pd.read_csv(measured_displacement_file,sep=r"\s+",skiprows=1,names=["days", "EXT-K1", "EXT-K2", "EXT-K3", "EXT-K4", "EXT-K5", "EXT-K6"])
measured_displacement["EXT-K6"] = (measured_displacement["EXT-K6"] - measured_displacement["EXT-K6"].iloc[0])
measured_displacement["time"] = pd.to_datetime(measured_displacement["days"])
measured_displacement["time"] = (measured_displacement["time"] - measured_displacement["time"].iloc[0]).dt.days

t_meas = measured_displacement["time"].to_numpy()
y_meas = measured_displacement["EXT-K6"].to_numpy()

# read simulation result OCR=1.5-----
displacement_ocr15_file = "output/ocr15/displacement.dat"
displacement_ocr15 = pd.read_csv(displacement_ocr15_file,sep=r"\s+",skiprows=1,names=["days", "displacement"])
displacement_ocr15["time"] = pd.to_datetime(displacement_ocr15["days"])
displacement_ocr15["time"] = (displacement_ocr15["time"] - displacement_ocr15["time"].iloc[0]).dt.days

# read simulation result OCR=1.7-----
displacement_ocr17_file = "output/ocr17/displacement.dat"
displacement_ocr17 = pd.read_csv(displacement_ocr17_file,sep=r"\s+",skiprows=1,names=["days", "displacement"])
displacement_ocr17["time"] = pd.to_datetime(displacement_ocr17["days"])
displacement_ocr17["time"] = (displacement_ocr17["time"] - displacement_ocr17["time"].iloc[0]).dt.days

# read simulation result OCR=2.0-----
displacement_ocr20_file = "output/ocr20/displacement.dat"
displacement_ocr20 = pd.read_csv(displacement_ocr20_file,sep=r"\s+",skiprows=1,names=["days", "displacement"])
displacement_ocr20["time"] = pd.to_datetime(displacement_ocr20["days"])
displacement_ocr20["time"] = (displacement_ocr20["time"] - displacement_ocr20["time"].iloc[0]).dt.days


t_sim = displacement_ocr15["time"].to_numpy()
y15_sim = displacement_ocr15["displacement"].to_numpy()
y17_sim = displacement_ocr17["displacement"].to_numpy()
y20_sim = displacement_ocr20["displacement"].to_numpy()

y15_sim_interp = np.interp(t_meas, t_sim, y15_sim)
mask = np.isfinite(y_meas) & np.isfinite(y15_sim_interp)
y_meas_valid = y_meas[mask]
y15_sim_valid = y15_sim_interp[mask]
r15, p15 = pearsonr(y_meas_valid, y15_sim_valid)

y17_sim_interp = np.interp(t_meas, t_sim, y17_sim)
mask = np.isfinite(y_meas) & np.isfinite(y17_sim_interp)
y_meas_valid = y_meas[mask]
y17_sim_valid = y17_sim_interp[mask]
r17, p17 = pearsonr(y_meas_valid, y17_sim_valid)

y20_sim_interp = np.interp(t_meas, t_sim, y20_sim)
mask = np.isfinite(y_meas) & np.isfinite(y20_sim_interp)
y_meas_valid = y_meas[mask]
y20_sim_valid = y20_sim_interp[mask]
r20, p20 = pearsonr(y_meas_valid, y20_sim_valid)


#-----------------------------------
#displacement
plt.figure()
plt.plot(measured_displacement["days"],measured_displacement["EXT-K6"],color="black",label="Extensometer measurement")

plt.plot(displacement_ocr15["days"],displacement_ocr15["displacement"],color="red",label="SubsApp OCR = 1.5")
plt.plot(displacement_ocr17["days"],displacement_ocr17["displacement"],color="lightskyblue",label="SubsApp OCR = 1.7")
plt.plot(displacement_ocr20["days"],displacement_ocr20["displacement"],color="green",label="SubsApp OCR = 2.0")

plt.title("Simulated Surface Vertical Displacement vs. Extensometer Measurements at the Nieuwolda Clay Site")
plt.xlabel(r"Time [day]")
plt.ylabel(r"Vetical displacement [mm]")
plt.xticks(rotation=45)
plt.legend()
textstr = (
    rf"$OCR = 1.5: R = {r15:.3f}$" "\n"
    rf"$OCR = 1.7: R = {r17:.3f}$" "\n"
    rf"$OCR = 2.0: R = {r20:.3f}$" "\n"
)
plt.text(0.05,0.05,textstr,transform=plt.gca().transAxes,fontsize=11,
    verticalalignment="bottom",bbox=dict(boxstyle="round",facecolor="white",alpha=0.8))

plt.grid(True,which="both",linestyle="--",alpha=0.5)
plt.gca().xaxis.set_major_locator(plt.MaxNLocator(8))
plt.tight_layout()
plt.savefig( "output/combined_displacement.png",dpi=300)
plt.close()


