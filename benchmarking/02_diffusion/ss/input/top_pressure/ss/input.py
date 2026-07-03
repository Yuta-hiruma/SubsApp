import pandas as pd

file_path = "top_pressure_swp.dat"  # Update with your actual file path

df = pd.read_csv(file_path, delim_whitespace=True, skiprows=1, names=["date", "pressure"])

# Convert to row format
bctop_time = df["date"].tolist()
bctop_pressure = df["pressure"].tolist()

print(bctop_time)
