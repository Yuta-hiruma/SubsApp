import pandas as pd
import matplotlib.pyplot as plt


# read measurement data-----
measure_displacement = "sub-cr_result/sub-cr_ch5.csv"
displacement_bp = pd.read_csv(measure_displacement, sep=",", skiprows=3, names=["time","ocr15","ocr17","ocr20"])

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
plt.plot(displacement_bp["time"],displacement_bp["ocr15"],marker="s",linestyle=":",color="darkred",label="Sub-CR OCR = 1.5")
plt.plot(displacement_bp["time"],displacement_bp["ocr17"],marker="s",linestyle=":",color="royalblue",label="Sub-CR OCR = 1.7")
plt.plot(displacement_bp["time"],displacement_bp["ocr20"],marker="s",linestyle=":",color="seagreen",label="Sub-CR OCR = 2.0")

plt.plot(displacement_ocr15["time"],displacement_ocr15["displacement"]/1000,color="red",label="SubsApp OCR = 1.5")
plt.plot(displacement_ocr17["time"],displacement_ocr17["displacement"]/1000,color="blue",label="SubsApp OCR = 1.7")
plt.plot(displacement_ocr20["time"],displacement_ocr20["displacement"]/1000,color="green",label="SubsApp OCR = 2.0")

plt.title("Comparison between Sub-CR and SubsApp")
plt.xlabel(r"Time $\psi$ [day]")
plt.ylabel(r"Vetical displacement [m]")
plt.legend()
plt.grid(True,which="both",linestyle="--",alpha=0.5)
plt.tight_layout()
plt.savefig( "output/combined_displacement.png",dpi=300)
plt.close()


