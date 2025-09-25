# SPyH FluidSim

![Python Version](https://img.shields.io/badge/python-3.9+-blue.svg)
![License](https://img.shields.io/badge/License-MIT-yellow.svg)

[cite_start]SPyH FluidSim is a lightweight, 2D fluid simulation package based on the Smoothed-Particle Hydrodynamics (SPH) method[cite: 1]. It's designed to be a straightforward tool for simulating and visualizing fluid dynamics. The project's full name is "SPyH FluidSim".



---

## Key Features

* **2D SPH Engine**: Implements the core SPH algorithm for fluid dynamics, calculating density and velocity changes for each particle.
* **Multiple Integration Schemes**: Includes Simple Euler (`simple1`), a predictor-corrector (`simple2`), and Leapfrog (`leapfrog`) integration methods to advance the simulation in time.
* **Selectable Kernel Functions**: Supports B-spline and Gaussian kernel functions for particle interaction calculations.
* [cite_start]**Customizable Physics**: Easily configure physical constants such as gravity, viscosity, and fluid density through the `PhysicsAttributes` class.
* **Efficient Neighbor Search**: Uses `scipy.spatial.KDTree` for fast and efficient particle neighbor finding within a given support radius.
* **Built-in Visualization**: Generates and saves simulation animations directly using `matplotlib`.

---

## Installation

To get started, clone the repository and install the required dependencies.

```bash
# Clone the repository
git clone [https://github.com/YoYo-XYZ/SPyH-FluidSim.git](https://github.com/YoYo-XYZ/SPyH-FluidSim.git)
cd SPyH-FluidSim

# Install dependencies (numpy, scipy, matplotlib are used)
pip install numpy scipy matplotlib