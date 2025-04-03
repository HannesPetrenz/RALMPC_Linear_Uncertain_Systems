# Robust Adaptive Learning Model Predictive Control (RALMPC) for Linear Systems with Parametric Uncertainties

## Overview
This repository contains MATLAB scripts for implementing **Robust Adaptive Learning Model Predictive Control (RALMPC)**. The approach is designed to handle **linear systems with parametric uncertainties**, combining **robust model predictive control (MPC)** with **adaptive learning techniques**. The primary focus is on improving control performance under uncertainty by leveraging iterative learning and parameter adaptation.

## Features
- **System Definition**: Defines system dynamics with **parametric uncertainties**.
- **Disturbance Modeling**: Configures **bounded disturbance sets**.
- **Cost Function Setup**: Defines quadratic cost functions for state tracking and control effort.
- **State & Input Constraints**: Implements polyhedral state and input constraints.
- **Adaptive Learning**: Initializes a **hypercube parameter set** for iterative parameter estimation.
- **Offline Computation**: Precomputes necessary parameters for **RALMPC** and computes an initial feasible solution.
- **Simulation for Different Horizons**: Runs **RALMPC and RLMPC** (Robust Learning MPC) for different **prediction horizons (N)**.
- **Result Storage**: Saves computed solutions for further analysis.

## Steps for Execution
1. **Perform Offline Computation**
   - Compute required parameters for RALMPC.
2. **Compute Initial Feasible Solution**
   - Generate an initial sample set.
3. **Solve Robust Optimal Control Problem**
   - Solve the robust optimal solution using the true system parameters.
4. **Run Robust Adaptive MPC (RAMPC)**
   - Initialize and solve **RAMPC** with adaptive parameter estimation.
5. **Run RALMPC and RLMPC for Different Horizons**
   - Iterate over different **prediction horizons** and compare results.
6. **Store Results**
   - Saves results in `.mat` files for further analysis.
## Folder Structure
This repository contains the following files:

- **`main.m`**: The main script that coordinates the execution.

- **`offlineComputation.m`**: Performs the offline computation of necessary parameters for RALMPC.

- **`get_Initalsolution.m`**: Generates an initial feasible solution for the RALMPC.

- **`solve_OS.m`**: Solves the robust optimal control problem using the true system parameters.

- **`solve_RAMPC.m`**: Solves the Robust Adaptive Model Predictive Control (RAMPC) problem with adaptive parameter estimation.

- **`solve_RALMPC.m`**: Solves the Robust Adaptive Learning Model Predictive Control (RALMPC) problem.

- **`generatedPlot.m`**: Generates plots for visualizing the results of the RALMPC and RAMPC simulations. This includes performance metrics and control trajectories for comparison.
  
## Installation & Requirements
- **MATLAB** (Tested on MATLAB R2023a and later)
- **MPT3 Toolbox** (for polyhedral operations) - [Homepage](https://www.mpt3.org/)
- **CasADi** (for optimization and nonlinear control) - [Homepage](https://web.casadi.org/)

## Usage
1. Clone the repository:
   ```bash
   git clone https://github.com/your-repo/RALMPC.git
   cd RALMPC

## License

MIT License

Copyright (c) 2025 Hannes Petrenz and University of California, Berkeley

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.