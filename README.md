# UAV Path Planning Benchmark

This repository contains the MATLAB source code associated with the paper:  
**A Multi-objective Benchmark for UAV Path Planning with Baseline Results**  
*Authors: Daison Darlan, Oladayo S. Ajani, Anand Paul, Rammohan Mallipeddi*

The paper introduces a flexible framework for generating synthetic environments (urban, suburban, and mountainous) and presents adapted versions of several multi-objective evolutionary algorithms (MOEAs) used in the comparative analysis for UAV path planning. This repository provides the MATLAB code for:

- **Environment Generation:** MATLAB scripts to generate realistic simulation environments including terrain modeling (using sinusoidal and Gaussian-based formulations), building placement, and no-fly zone integration.
- **Algorithm Implementations:** Adapted and benchmarked MATLAB implementations of several MOEAs for UAV path planning, including:
  - NSGA-II and NSGA-II/SDR
  - MOEA/D and MOEA/D-AWA
  - HypE
  - 𝐼𝑆𝐷𝐸+
  - MOEA-2DE

## Getting Started

### Prerequisites
- MATLAB R2019b or later
- Required MATLAB toolboxes

### Installation
1. Clone the repository:
   ```bash
   git clone https://github.com/your-username/UAV-Path-Planning-Benchmark.git
2. Open the cloned repository in MATLAB.

### Running the Code
- Generate Environments: In MATLAB, run the environment generator script to create and visualize sample environments:
  ```bash
  run('examples/run_environment.m')
- Run Benchmark Experiments: Execute the benchmark script to run comparative experiments on the provided test scenarios:
  ```bash
  run('examples/run_benchmark.m')
