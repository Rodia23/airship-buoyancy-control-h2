# Regenerative Hydrogen Buoyancy Control for Heavy-Lift Airships

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[INSERT ZENODO DOI BADGE HERE]

This repository contains the simulation source code, datasets, and LaTeX manuscript for the paper: **"Regenerative Hydrogen Buoyancy Control for Heavy-Lift Airships: Actuator Requirements and Scaling Effects"**.

The project models a closed-loop hydrogen regenerative cycle for buoyancy control during heavy-lift payload release operations. It uses a Proportional-Derivative (PD) controller coupled with flight dynamics and thermodynamic models to determine the minimum actuator capacity required to stabilize representative airships (e.g., Pathfinder 3, Flying Whales LCA60T, Zeppelin NT, and LZ 129 Hindenburg).

## Repository Structure

```text
airship-buoyancy-h2/
│
├── paper_manuscript/         # LaTeX source files and figures for the MDPI submission
│   ├── articulo_mdpi.tex
│   └── *.png
│
├── simulation_code/          # Python source code for dynamic simulation
│   ├── parametros.py         # System constants, airship specs, and controller gains
│   ├── simulacion.py         # Base ODE solver for flight dynamics
│   ├── t2_simulacion.py      # Base mission and multi-drop stress test simulation
│   ├── t3_barrido.py         # Parametric sweep for minimum actuator requirements
│   ├── ciclo_balance.py      # Energy balance, condensation time, and solar regeneration
│   ├── comparacion_O2.py     # Mass balance comparison (onboard O2 vs. atmospheric O2)
│   ├── generar_figura9.py    # Formatting script for condensation power plot
│   └── analisis_optimos.py   # Data extraction tool for minimum and optimal rates
│
├── LICENSE
└── README.md
```
## Requirements
To run the simulations and reproduce the plots from the paper, you need Python 3.x and the following libraries:

## Bash
pip install numpy scipy matplotlib
Usage and Reproduction
All scripts are configured to run out-of-the-box and generate the corresponding .csv datasets and high-resolution .png figures used in the manuscript.

Configure Parameters: Modify parametros.py if you wish to test a custom airship configuration.

Run the Parametric Sweep: ```bash
python simulation_code/t3_barrido.py

*This evaluates the minimum H2 combustion rates across the 4 representative airships.*
Run the Dynamic Mission (Stress Test):

Bash
python simulation_code/t2_simulacion.py
This simulates the transient stabilization during a multi-drop payload release.

Run the Energy Balance:

Bash
python simulation_code/ciclo_balance.py
Calculates the solar regeneration capabilities and daily operational autonomy.

## Citation
If you use this code or our methodology in your research, please cite our work:

Espinoza Leon, A.E.; Ortiz Porras, J.E. (2026). Regenerative Hydrogen Buoyancy Control for Heavy-Lift Airships: Actuator Requirements and Scaling Effects. MDPI [Journal Name - PENDING].

## License
This project is licensed under the MIT License - see the LICENSE file for details.