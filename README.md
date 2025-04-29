---

# 🧪 SPyH Fluid Simulation (Numpy)

This project is a 2D fluid dynamics simulator based on **Smoothed Particle Hydrodynamics (SPH)**. It features a real-time interactive interface using **Tkinter** and visualizes particle motion with **Matplotlib**.

---

## ✨ Features

- Customizable particle setup (position, count, box size)
- Adjustable physical parameters: viscosity, gravity, smoothing radius
- Realistic fluid behavior using SPH principles
- Collision handling with bounding box walls
- Live rendering and setup preview
- Linked-list based neighbor search optimization

---

## 🖼️ GUI Preview

---

## 📦 Project Structure

```
.
├── main.py                 # Entry point: launches the GUI
├── gui.py                  # Tkinter-based interface for parameter control
├── simulation_backend.py  # Core SPH logic and physics engine
├── illustration.py         # (Expected to contain Matplotlib plotting functions)
```

---

## 🛠️ Installation

1. **Clone the repo**:

```bash
git clone https://github.com/yourusername/sph-fluid-simulation.git
cd sph-fluid-simulation
```

2. **Install dependencies**:

```bash
pip install numpy matplotlib
```

---

## 🚀 Run the Simulation

```bash
python main.py
```

The GUI will launch. You can then:

- Set your simulation box size
- Choose initial particle amount and distribution
- Adjust smoothing radius, viscosity, gravity, frame count, and camera speed
- Click **"Show Setup"** to preview
- Click **"Perform Simulation"** to run

---

## 🧠 Simulation Method

This simulator uses:

- SPH equations with B-spline kernel function
- Linked-list cell-based neighbor search
- Euler integration with adjustable substeps
- Symmetrized forces for realistic pressure and viscosity dynamics
- 2D rectangular boundary with collision damping

---

## 📃 License

MIT License.

---