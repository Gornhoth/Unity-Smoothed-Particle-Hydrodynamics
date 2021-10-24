# Unity-Smoothed-Particle-Hydrodynamics
This repository contains two Unity projects showcasing fluid simulations with the SPH method. There are three different implementations using:
1. MonoBehaviour
2. ComputeShader
3. Entity Component System (DOTS)

As spatial acceleration structure a uniform grid (as described by Simon Green) has been chosen. Particles are represented as uniformely sized spheres. There is no collision geometry involved (only boundary condition for simulation space; particles are guided only by SPH). The implementations are deliberately similar to ease understanding of the SPH method and its performance impact when implemented using different parallelisation techniques, namely the Entity-Component-System (ECS) and ComputeShaders.

(placeholder for a bit more detailed description of the two unity projects, especially the one containing MonoBehaviour and ComputeShader as the mono behaviour one is deactivated by default)

Knowledge and inspiration for the code in this repository was taken from:
- Doyub Kim https://github.com/doyubkim/fluid-engine-dev
- Müller, Matthias, David Charypar, and Markus Gross. "Particle-based fluid simulation for interactive applications."
- Simon Green "Particle simulation using cuda."
- https://lucasschuermann.com/writing/implementing-sph-in-2d
- Marcel Kießlich
- Pascal Blessing

(placeholder for future images/videos of the simulation)
