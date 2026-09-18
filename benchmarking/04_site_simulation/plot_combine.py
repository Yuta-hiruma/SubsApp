import pandas as pd
import matplotlib.pyplot as plt


#read measurement -----------
measured_displacement_file = "measurement/lvl.dat"
measured_displacement = pd.read_csv(measured_displacement_file,sep=r"\s+",skiprows=1,names=["days", "EXT-K1", "EXT-K2", "EXT-K3", "EXT-K4", "EXT-K5", "EXT-K6"])
measured_displacement["EXT-K6"] = (measured_displacement["EXT-K6"] - measured_displacement["EXT-K6"].iloc[0])
measured_displacement["time"] = pd.to_datetime(measured_displacement["days"])
measured_displacement["time"] = (measured_displacement["time"] - measured_displacement["time"].iloc[0]).dt.days


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

#-----------------------------------
#displacement
plt.figure()
plt.plot(measured_displacement["time"],measured_displacement["EXT-K6"],color="black",label="Extensometer measurement")

plt.plot(displacement_ocr15["time"],displacement_ocr15["displacement"],color="red",label="SubsApp OCR = 1.5")
plt.plot(displacement_ocr17["time"],displacement_ocr17["displacement"],color="blue",label="SubsApp OCR = 1.7")
plt.plot(displacement_ocr20["time"],displacement_ocr20["displacement"],color="green",label="SubsApp OCR = 2.0")

plt.title("Simulated Surface Vertical Displacement vs. Extensometer Measurements at the Nieuwolda Clay Site")
plt.xlabel(r"Time [day]")
plt.ylabel(r"Vetical displacement [mm]")
plt.legend()
plt.grid(True,which="both",linestyle="--",alpha=0.5)
plt.gca().xaxis.set_major_locator(plt.MaxNLocator(8))
plt.tight_layout()
plt.savefig( "output/combined_displacement.png",dpi=300)
plt.close()


