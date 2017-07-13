# Active-Brownian-
Active Brownian particle movement and their cluster formation

MD simulation of self propulling Brownian particles 

Number of particles used = 4096
dt = 0.00001
Time steps for equilibriation = 1000000
Total time steps = 2000000
Density of a single ensemble = 0.80
vp : Propulsion Velocity

The code is parallelized using OpenMP

Brownian particles interact with each other via the modified Lennard-Jones potential
known as the WCA potential
This is also called soft active Browian simulation
Cut off distance for the potential is 2^(1/6)


Trajectory of the particle is determined using Langevin equations
