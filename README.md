# Unity Smoothed Particle Hydrodynamics
This repository contains two Unity projects showcasing fluid simulations with the SPH method. There are three different implementations using:
1. MonoBehaviour (single-threaded CPU)
2. ComputeShader (GPU)
3. Entity Component System (DOTS, multi-threaded CPU via JobSystem)

As spatial acceleration structure a uniform grid (as described by Simon Green) has been chosen. Particles are represented as uniformely sized spheres. Their size (= interaction radius) has been chosen with respect to the grid cell size (= 2x particle radius) so that one particle must only consider 8 cells in total for its neighbour search. There is no collision geometry involved (only boundary condition for simulation space; particles are guided only by SPH). The implementations are deliberately similar to ease understanding of the SPH method and its performance impact when implemented using different parallelisation techniques, namely the Entity-Component-System (ECS) and ComputeShaders.

Usage: Open one of the two unity projects. After restoring packages you may start the application in the editor (you can also build although note the difference of building for ECS via BuildConfiguration). The MonoBehaviour implementation is disabled by default as it shares a scene with the ComputeShader implementation that is enabled by default (see "MonoBehaviour CPU SPH" and "ComputeShader GPU SPH" gameobjects in "SPH" scene). You can influence simulation parameters and view debug information on the gameobjects' SPH management scripts (ParticleManager.cs for MonoBehaviour implementation and ComputeShaderParticleManager.cs for ComputeShader implementation). Be careful with the MonoBehaviour implementaion's simulation parameters, as you can easily freeze your unity editor if you want to simulate too many particles. The ECS project does not have any way to change and view simulation parameters in the editor. Therefore, you can only influence the simulation by editing the SPHSystem.cs file. Also be careful when changing the particle interaction radius and grid cell size because it could break the neighbour search.

Known limitations:
- Due to the hashgrid's limitation of having a maximum number of particles per cell, the simulation may become unstable if more particles than that occupy a cell's space. See issue #1 for more information.

Knowledge and inspiration for the code in this repository was taken from:
- Doyub Kim https://github.com/doyubkim/fluid-engine-dev
- Müller, Matthias, David Charypar, and Markus Gross. "Particle-based fluid simulation for interactive applications."
- Simon Green "Particle simulation using cuda."
- https://lucasschuermann.com/writing/implementing-sph-in-2d
- Marcel Kießlich
- Pascal Blessing

https://user-images.githubusercontent.com/34295773/138610298-1851bf02-e4a4-4335-a02d-151825802617.mp4

![sensibleConfig](https://user-images.githubusercontent.com/34295773/138756364-eb689db3-0481-45ed-b74c-f559fe62911a.png)

FYI: If you ever wondered just how big the difference between (parallel-)CPU and GPU implementations is, here is a plot measuring the average frametimes over the number of particles for each implementation on a system with an AMD Ryzen 3600 and NVIDIA GeForce RTX 2060 Super:

![AvgFrametimeBase_System2](https://user-images.githubusercontent.com/34295773/138755834-7dfae6c7-9682-43d3-b97d-67c051705f5b.png)
