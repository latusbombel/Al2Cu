# Al2Cu Structure Modeling & Energy Minimization

> ⚠️ **Project Status: Active Development**  
> This repository is currently under active development. Core simulation routines, minimization algorithms, and analysis scripts are being refined. APIs, file formats, and CLI interfaces may change without notice.

---

## Overview

This project provides computational tools for modeling, relaxing, and analyzing intermetallic **$\text{Al}_2\text{Cu}$** ($\theta$-phase) crystal structures. The codebase combines C++ for numerical algorithms (geometry optimization, energy evaluation, and X-ray diffraction calculation) and Python for configuration generation and post-processing visualizations.

---

## Repository Structure

```text
├── .vscode/
│   └── settings.json           # VS Code workspace settings
├── Al2Cu.cif                   # Crystallographic Information File for Al2Cu
├── classes.h                   # Class declarations (lattice, atoms, simulation box)
├── functions.cpp               # Supporting mathematical and helper routines
├── Al2Cu.cpp                   # Main program entry point for simulation routines
├── minimalizacja.cpp           # Energy minimization and geometry relaxation routines
├── xrd.cpp                     # C++ engine for X-ray diffraction (XRD) spectrum calculations
├── poz_startowe.py             # Python script for generating initial atom positions
├── plot_minimalization.py      # Plotting convergence and energy curves
├── xrd.py                      # Python visualization script for XRD spectra
├── start_pos.dat               # Coordinate data file containing initial atomic positions
└── wyniki/                     # Simulation results and output artifacts
    ├── dane.dat                # Numerical summary data
    ├── minimalizacja.dat       # Energy vs. iteration step data
    ├── energia2.png            # Plot of potential energy minimization
    ├── polozenie_pocz.dat/.xyz # Initial atomic coordinates (DAT and XYZ formats)
    ├── polozenie_min.dat/.xyz  # Coordinates at energy minimum
    ├── polozenie_end.dat       # Final coordinates after simulation run
    ├── xrd.dat                 # Computed 2-theta vs. intensity XRD data
    ├── xrd_pocz.png            # XRD pattern of the starting configuration
    └── xrd_end.png             # XRD pattern of the relaxed/final configuration
```


---

## Features

1. **Lattice & Structure Modeling**:
   - Reference crystallographic parameters defined in `Al2Cu.cif`.
   - Starting configuration generation with periodic boundary conditions via `poz_startowe.py`.

2. **Energy Minimization**:
   - Iterative relaxation of atomic coordinates implemented in `minimalizacja.cpp`.
   - Convergence tracking and trajectory output (`polozenie_pocz.xyz` $\rightarrow$ `polozenie_min.xyz`).
   - Automated plotting of the minimization trajectory with `plot_minimalization.py`.

3. **Structural Characterization (XRD)**:
   - Calculation of kinematic powder diffraction patterns using atomic scattering factors in `xrd.cpp`.
   - Comparative spectrum generation (`xrd_pocz.png` vs. `xrd_end.png`) through `xrd.py` to verify phase stability and structural relaxation.

---

## Prerequisites

- **C++ Compiler**: GCC ($\ge 9.0$) or Clang ($\ge 10.0$) supporting C++17.
- **Python**: Version 3.8+ with the following packages:
  ```bash
  pip install numpy matplotlib scipy
  ```
- **Structure Visualization** *(Optional)*: [OVITO](https://www.ovito.org/) or [VESTA](https://jp-minerals.org/vesta/en/) to inspect `.xyz` and `.cif` trajectory files.

---

## Compilation & Usage

### 1. Build C++ Executables

Compile the source routines using standard optimization flags:

```bash
# Compile energy minimization routine
g++ -O3 -std=c++17 minimalizacja.cpp functions.cpp -o minimalizacja

# Compile XRD calculation routine
g++ -O3 -std=c++17 xrd.cpp functions.cpp -o xrd
```

### 2. Generate Initial Positions

Generate the starting coordinates if rebuilding the domain:

```bash
python poz_startowe.py
```

### 3. Run Energy Minimization

Execute the relaxation binary to produce output coordinate and energy files in `wyniki/`:

```bash
./minimalizacja
```

Plot the resulting energy progression:

```bash
python plot_minimalization.py
```

### 4. Compute and Plot XRD Patterns

Run the diffraction simulation and plot the comparative spectra:

```bash
./xrd
python xrd.py
```

---

## Output Data Format

- **`wyniki/minimalizacja.dat`**: Space-separated values tracking iteration step, total potential energy, and maximum interatomic force component.
- **`wyniki/polozenie_*.xyz`**: Standard XYZ trajectory format containing chemical species (`Al`, `Cu`) and Cartesian coordinates in Angstroms ($x$, $y$, $z$).
- **`wyniki/xrd.dat`**: Two-column array containing scattering angle $2\theta$ (degrees) and relative intensity $I$ (arbitrary units).

---

## Roadmap / Planned Improvements

- [ ] Add support for temperature-dependent molecular dynamics (NVT/NPT ensembles).
- [ ] Include automated unit and regression tests for energy minimization routines.
- [ ] Add MSD calculation.
