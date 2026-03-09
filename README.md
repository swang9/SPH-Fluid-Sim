# SPH 2D Fluid Simulation

A real-time 2D fluid simulation built in C++ using Smoothed Particle Hydrodynamics (SPH), with OpenGL for rendering and Eigen for linear algebra operations.

https://github.com/user-attachments/assets/39578cd2-44ee-49c7-a032-425618f7540c


## Overview

SPH is a particle-based method for simulating fluid dynamics. Each particle carries physical properties like position, velocity, density, and pressure. At each timestep, forces are computed between neighboring particles using smoothing kernels, producing realistic fluid-like behavior.

This project implements the core SPH algorithm from scratch, including:
- Density and pressure estimation per particle
- Pressure and viscosity force calculations
- Boundary collision handling
- Real-time OpenGL visualization

## Tech Stack

- **Language:** C++
- **Rendering:** OpenGL, GLUT
- **Linear Algebra:** Eigen
- **Build System:** Makefile with OS-detection for macOS and Linux

## Getting Started

### Prerequisites

Install OpenGL and GLUT:
```
sudo apt-get install freeglut3-dev
```

Install Eigen:
```
sudo apt-get install libeigen3-dev
```

### Building and Running

1. Clone the repository
   ```
   git clone https://github.com/swang9/SPH-Fluid-Sim.git
   cd SPH-Fluid-Sim
   ```

2. In the Makefile, comment and uncomment the section matching your OS (macOS or Linux)

3. If on Linux, change `arc4random` to `rand` in `main.cpp`

4. Build and run
   ```
   make
   ./sph
   ```

> Note: currently supported on macOS and Linux only.

<!--
## Instructions

For MacOS and Linux.

Install OpenGL and GLUT libraries with

`sudo apt-get install freeglut3-dev`

Install Eigen library with

`sudo apt-get install libeigen3-dev`

Based on your operating system, comment and uncomment the relevant sections in the Makefile.

If on Linux, change `arc4random` in `main.cpp` to `rand`.


`make`, then run with `./sph`.
-->
