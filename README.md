# SPH 2D Fluid Simulation

## Instructions

For MacOS and Linux.

Install OpenGL and GLUT libraries with

`sudo apt-get install freeglut3-dev`

Install Eigen library with

`sudo apt-get install libeigen3-dev`

Based on your operating system, comment and uncomment the relevant sections in the Makefile.

If on Linux, change `arc4random` in `main.cpp` to `rand`.

`make`, then run with `./sph`.