# SubsApp

**SubsApp** is a one-dimensional hydro-mechanical model developed to simulate shallow land subsidence driven by pore-water pressure changes in saturated and unsaturated soils. The model incorporates the **Soft Soil Creep (SSC)** formulation to simulate time-dependent vertical soil deformation and accounts for matric suction in the unsaturated zone. SubsApp is designed for investigating shallow subsurface deformation in soft soils, particularly peat and clay.

SubsApp has been tested through several validation and benchmarking applications. However, errors or limitations not identified during these tests may emerge in future applications. Users are encouraged to report any errors encountered in the model or its documentation.

The **user manual**, available in the `doc/` directory, provides detailed information on the model formulation, numerical implementation, and usage.

## Features

- **One-dimensional hydro-mechanical coupling** between pore-water pressure and vertical soil deformation.
- **Layered soil profiles** to represent different soil types and stratigraphic conditions.
- **Time-dependent soil deformation** based on the Soft Soil Creep (SSC) model.
- **Pore-water pressure simulation** using a one-dimensional diffusion equation, including matric suction in the unsaturated zone.
- **Variable hydraulic conductivity** in unsaturated soils based on the Mualem–van Genuchten formulation.
- **Flexible top boundary conditions**, allowing either prescribed pore-water pressure or water flux.
- **Vertical displacement simulation** to estimate soil-column deformation and land subsidence.

## Project Structure

```text
SubsApp/
├── data/             # Input data and soil material parameters
├── doc/              # User manual and documentation
├── results/          # Simulation outputs
├── src/              # Model source code
│   └── input/        # Initial and boundary conditions
├── benchmarking/     # Benchmarking tests and associated input data
└── README.md
```

## Software Requirements

SubsApp has been tested using **Python 3.11.5**.

Python can be downloaded from the [official Python website](https://www.python.org/).

The following software and package versions were used in the tested environment:

| Software / Package | Version |
|---|---|
| Python | 3.11.5 |
| geopandas | 0.14.4 |
| matplotlib | 3.7.2 |
| matplotlib-inline | 0.1.6 |
| numba | 0.61.2 |
| numba_celltree | 0.4.1 |
| numpy | 1.26.4 |
| numpydoc | 1.5.0 |
| pandas | 2.3.1 |
| scipy | 1.11.1 |
| GMT | 6.6.0 |

**Note:** These versions represent the tested software environment. Compatibility with other versions has not been verified.

## Usage

### 1. Prepare Input Data

The model input files are located in the `data/` directory.

#### `input.py`

This file defines the simulation setup, including:

- Site name
- Simulation period
- Initial groundwater level
- Soil stratigraphy
- Soil type assigned to each layer
- Boundary condition types

#### `input_soil_material_parameters.py`

This file defines the material properties for each soil type, including:

- Soil compression parameters
- Hydraulic parameters
- Soil unit weights

### 2. Prepare Initial and Boundary Conditions

The initial and boundary condition files are located in the `src/input/` directory.

```text
src/
└── input/            # Initial and boundary condition files
```

These files provide the initial and boundary conditions required for the simulation.

For detailed instructions on preparing input data and implementing the model, refer to the user manual in the `doc/` directory.

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

## Documentation

Further details on the theoretical formulation, numerical implementation, model parameters, and validation tests are provided in the user manual.

Users are encouraged to consult the documentation before applying SubsApp to new study areas or soil conditions.

## Reporting Issues

Although SubsApp has undergone validation and benchmarking, errors or limitations may remain. Users are encouraged to report bugs, unexpected model behaviour, or documentation errors through the GitHub repository's **Issues** section.


## Author

Yuta Hiruma

## License

This software is developed at Utrecht University. 
