# SubsApp

**SubsApp** is a one-dimensional hydro-mechanical model for simulating shallow. The model is based on the **Soft Soil Creep (SSC)** fomulation and includes matric suction in unsaturated zone. It is designed for shallow subsurface deformation analysis of soft soils such as peat and clay.

## Features

* One-dimensional hydro-mechanical coupling
* Layered soil profile representation
* Subsidence simulation using the Soft Soil Creep (SSC) model
* Simulate pore pressure using defusion equation.
* Hydraulic conductivity in unsaturated zone is updated by Mualem-van Genuchten formulation
* Matric suction is calcurated
* Pore pressure boundary conditions is to be derived from SWAP simulations
* Surface displacement output

## Project Structure

```text
SubsApp/
├── data/           # Input files and soil parameters
├── doc/            # Documentation
├── results/        # Simulation results
├── src/            # Script code
└── README.md
```

## Usage

### 1. Prepare Input Data

Input files are located in the `data/` directory.

#### `input.py`

Define:

* Site name
* Simulation period
* Initial groundwater level
* Soil stratigraphy
* Soil type assigned to each layer

#### `input_soil_material_parameters.py`

Define the material properties for each soil type, including:

* Compression parameters
* Hydraulic parameters
* Creep parameters
* Unit weights

### 2. Generate Boundary Conditions

Top and bottom pore pressure boundary conditions can be generated using the SWAP model.

Run SWAP:

```bash
cd swap
./swap420
```

Process SWAP outputs:

```bash
cd output_swap

python3 pressure_bc.py
python3 pressure_ini.py
python3 water_content.py
```
This converts the swap output to SubsApp's boundary condition.

### 3. Run the Subsidence Model

```bash
cd src
python3 main.py
```

## Results

Simulation outputs are written to the `results/` directory and may include:

* Surface settlement time series
* Pressure distribution
* Effective saturation


## References



## Author

Yuta Hiruma

## License

This software is developed at Utrecht University. 
