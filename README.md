# SPyH FluidSim

![Python Version](https://img.shields.io/badge/python-3.9+-blue.svg)
![License](https://img.shields.io/badge/License-MIT-yellow.svg)

SPyH FluidSim is a lightweight, 2D fluid simulation package based on the Smoothed-Particle Hydrodynamics (SPH) method. It's designed to be a straightforward tool for simulating and visualizing fluid dynamics. The project's full name is "SPyH FluidSim".



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
```

## Quick Start
You can set up and run a simple "water droplet" simulation with just a few lines of code. The following example initializes particles in a rectangle, configures the simulation, and generates a GIF of the result.

from SPyH_FluidSim import Simulation, SimAttributes, ParticleAttributes, PhysicsAttributes, ParticleInitialize

```python

# 1. Define simulation attributes (boundary, timestep, duration)
sim_config = SimAttributes(bound=[1.0, 1.0], collision_coefficient=0.8, Dt=0.01, substep=5, totaltime=4.0)

# 2. Define physical properties of the fluid
physics_config = PhysicsAttributes(GRAVITY=9.81, VISCOSITY=0.1, WATER_DENSITY=1000.0)

# 3. Define particle and kernel properties
# The support radius should be larger than the particle gap.
particle_config = ParticleAttributes(support_radius=0.05, kernel_type='bspline')

# 4. Initialize particle positions and velocities
particles = ParticleInitialize(gap=0.02)
particles.rectangle(width=0.8, length=0.4, center=[-0.5, 0.0], velocity=[0.0, 0.0])

# 5. Initialize and run the simulation
particle_config._initialize_particle(particles, physics_config)
sim = Simulation(sim_config, particle_config)
sim.update_fulltime(step_method='leapfrog')

# 6. Create and save the animation
anim = sim.create_anim()
sim.save_anim(anim, "dam_break_simulation")
```
This script will produce an animation file named dam_break_simulation.gif in your project directory.

## License
This project is licensed under the MIT License. The software is provided "as is" without any warranty. In no event shall the authors be liable for any claim or damages.


Copyright (c) 2025 Thammatorn Jamraschai.