using System.IO;
using System.Runtime.InteropServices;
using Unity.Burst;
using Unity.Collections;
using Unity.Entities;
using Unity.Jobs;
using Unity.Rendering;
using UnityEngine;

// ReSharper disable InconsistentNaming
public class SPHSystem : SystemBase
{
    // ReSharper disable MemberCanBePrivate.Global
    public float radius = 1f; // Particle radius, cell size is 2 * particle Radius

    public int dimensions = 100;
    public int numberOfParticles = 10000;
    // Optimisation rule of thumb: (numberOfParticles^(1/3)) * 0.8
    public int maximumParticlesPerCell = 18;
    public float mass = 1f;
    public float gasConstant = 2000.0f;
    public float restDensity = 2f;
    public float viscosityCoefficient = 0.25f; 
    public Vector3 gravity = new Vector3(0.0f, -9.81f, 0.0f);
    public float damping = -0.5f;
    public float dt = 0.0008f;

    private NativeArray<Particle> _particles;
    public NativeMultiHashMap<int, int> _hashGrid;
    public NativeArray<int> _neighbourTracker;
    public NativeArray<int> _neighbourList; // Stores all neighbours of a particle aligned at 'particleIndex * maximumParticlesPerCell * 8'
    public NativeArray<float> densities;
    public NativeArray<float> pressures;
    public NativeArray<Vector3> forces;
    public NativeArray<Vector3> velocities;
    
    public Mesh particleMesh;
    public float particleRenderSize = 40f;
    public Material material;
    public float radius2;
    public float radius3;
    public float radius4;
    public float radius5;
    public float mass2;
    
    private ComputeBuffer _particleColorPositionBuffer;
    private ComputeBuffer _argsBuffer;
    private static readonly int SizeProperty = Shader.PropertyToID("_size");
    private static readonly int ParticlesBufferProperty = Shader.PropertyToID("_particlesBuffer");
    
    [Tooltip("The absolute accumulated simulation steps")]
    public int elapsedSimulationSteps;

    private int batchCount = 256;
    
    // ReSharper restore MemberCanBePrivate.Global
    // ReSharper restore InconsistentNaming

    [StructLayout(LayoutKind.Sequential, Size=28)]
    private struct Particle
    {
        public Vector3 Position;
        public Vector4 Color;
    }
    
    protected override void OnStartRunning()
    {
        Application.targetFrameRate = -1;
        var stuff = EntityManager.GetSharedComponentData<RenderMesh>(GetSingletonEntity<ParticlePrefabTag>());
        particleMesh = stuff.mesh;
        material = stuff.material;
        radius2 = radius * radius;
        radius3 = radius2 * radius;
        radius4 = radius3 * radius;
        radius5 = radius4 * radius;
        mass2 = mass * mass;
        
        RespawnParticles();
        InitNeighbourHashing();
        InitComputeBuffers();
    }

    protected override void OnStopRunning()
    {
        EntityManager.CompleteAllJobs();

        ReleaseNative();
    }

    private void ReleaseNative()
    {
        _particles.Dispose();
        _hashGrid.Dispose();
        _neighbourList.Dispose();
        _neighbourTracker.Dispose();
        densities.Dispose();
        pressures.Dispose();
        forces.Dispose();
        velocities.Dispose();

        _particleColorPositionBuffer.Dispose();
        _argsBuffer.Dispose();
    }

    #region Initialisation
    
    private void RespawnParticles()
    {
        _particles = new NativeArray<Particle>(numberOfParticles, Allocator.Persistent);
        densities = new NativeArray<float>(numberOfParticles, Allocator.Persistent);
        pressures = new NativeArray<float>(numberOfParticles, Allocator.Persistent);
        forces = new NativeArray<Vector3>(numberOfParticles, Allocator.Persistent);
        velocities = new NativeArray<Vector3>(numberOfParticles, Allocator.Persistent);

        int particlesPerDimension = Mathf.CeilToInt(Mathf.Pow(numberOfParticles, 1f / 3f));

        int counter = 0;
        while (counter < numberOfParticles)
        {
            for (int x = 0; x < particlesPerDimension; x++)
            for (int y = 0; y < particlesPerDimension; y++)
            for (int z = 0; z < particlesPerDimension; z++)
            {
                Vector3 startPos = new Vector3(dimensions - 1, dimensions - 1, dimensions - 1) - new Vector3(x / 2f, y / 2f, z / 2f)  /*- new Vector3(Random.Range(0f, 0.01f), Random.Range(0f, 0.01f), Random.Range(0f, 0.01f))*/;
                _particles[counter] = new Particle
                {
                    Position = startPos,
                    Color = Color.white
                };
                densities[counter] = -1f;
                pressures[counter] = 0.0f;
                forces[counter] = Vector3.zero;
                velocities[counter] = Vector3.zero;

                if (++counter == numberOfParticles)
                {
                    return;
                }
            }
        }
    }

    private void InitNeighbourHashing()
    {
        _neighbourList = new NativeArray<int>(numberOfParticles * maximumParticlesPerCell * 8, Allocator.Persistent); // 8 because we consider 8 cells
        _neighbourTracker = new NativeArray<int>(numberOfParticles, Allocator.Persistent);
        SpatialHashing.CellSize = radius * 2; // Setting cell-size h to particle diameter.
        SpatialHashing.Dimensions = dimensions;
        _hashGrid = new NativeMultiHashMap<int, int>(numberOfParticles, Allocator.Persistent);
    }

    void InitComputeBuffers()
    {
        uint[] args = {
            particleMesh.GetIndexCount(0),
            (uint) numberOfParticles,
            particleMesh.GetIndexStart(0),
            particleMesh.GetBaseVertex(0),
            0
        };
        _argsBuffer = new ComputeBuffer(1, args.Length * sizeof(uint), ComputeBufferType.IndirectArguments);
        _argsBuffer.SetData(args);
        
        _particleColorPositionBuffer = new ComputeBuffer(numberOfParticles, sizeof(float) * ( 3 + 4 ));
        _particleColorPositionBuffer.SetData(_particles);
    }
    
    #endregion
    
    [BurstCompile]
    private struct RecalculateHashGrid : IJobParallelFor
    {
        public NativeMultiHashMap<int, int>.ParallelWriter hashGrid; // Hash of cell to particle indices.
        [ReadOnly] public NativeArray<Particle> particles;

        [ReadOnly] public int dimensions;
        [ReadOnly] public float cellSize;

        public void Execute(int index)
        {
            hashGrid.Add(SpatialHashing.Hash(SpatialHashing.GetCell(particles[index].Position, cellSize), dimensions), index);
        }
    }
    
    [BurstCompile]
    private struct BuildNeighbourList : IJobParallelFor
    {
        [NativeDisableParallelForRestriction] public NativeArray<int> neighbourTracker;
        [NativeDisableParallelForRestriction] public NativeArray<int> neighbourList; // Stores all neighbours of a particle aligned at 'particleIndex * maximumParticlesPerCell * 8'
        [ReadOnly] public NativeMultiHashMap<int, int> hashGrid; // Hash of cell to particle indices.
        [ReadOnly] public NativeArray<Particle> particles;
        
        [ReadOnly] public float cellSize;
        [ReadOnly] public float radiusSquared;
        [ReadOnly] public int maximumParticlesPerCell;
        [ReadOnly] public int dimensions;

        public void Execute(int index)
        {
            neighbourTracker[index] = 0;
            var cell = SpatialHashing.GetCell(particles[index].Position, cellSize);
            var cells = GetNearbyKeys(cell, particles[index].Position);

            for (int j = 0; j < cells.Length; j++)
            {
                if (!hashGrid.ContainsKey(cells[j])) continue;
                var neighbourCell = hashGrid.GetValuesForKey(cells[j]);
                int counter = 0;
                foreach (var potentialNeighbour in neighbourCell)
                {
                    if (potentialNeighbour == index) continue;
                    if (( particles[potentialNeighbour].Position - particles[index].Position ).sqrMagnitude < radiusSquared) // Using squared length instead of magnitude for performance
                    {
                        neighbourList[index * maximumParticlesPerCell * 8 + neighbourTracker[index]++] = potentialNeighbour;
                        if (++counter == maximumParticlesPerCell) break;    // Prevent potential UB in neighbourList by only allowing maximumParticlesPerCell neighbours from one cell.
                    }
                }
            }
        }
        
        private NativeArray<int> GetNearbyKeys(Vector3Int originIndex, Vector3 position)
        {
            NativeArray<int> nearbyBucketIndicesX = new NativeArray<int>(8, Allocator.Temp);
            NativeArray<int> nearbyBucketIndicesY = new NativeArray<int>(8, Allocator.Temp);
            NativeArray<int> nearbyBucketIndicesZ = new NativeArray<int>(8, Allocator.Temp);
            for (int i = 0; i < 8; i++)
            {
                nearbyBucketIndicesX[i] = originIndex.x;
                nearbyBucketIndicesY[i] = originIndex.y;
                nearbyBucketIndicesZ[i] = originIndex.z;
            }

            if (( originIndex.x + 0.5f ) * cellSize <= position.x)
            {
                nearbyBucketIndicesX[4] += 1;
                nearbyBucketIndicesX[5] += 1;
                nearbyBucketIndicesX[6] += 1;
                nearbyBucketIndicesX[7] += 1;
            }
            else
            {
                nearbyBucketIndicesX[4] -= 1;
                nearbyBucketIndicesX[5] -= 1;
                nearbyBucketIndicesX[6] -= 1;
                nearbyBucketIndicesX[7] -= 1;
            }

            if (( originIndex.y + 0.5f ) * cellSize <= position.y)
            {
                nearbyBucketIndicesY[2] += 1;
                nearbyBucketIndicesY[3] += 1;
                nearbyBucketIndicesY[6] += 1;
                nearbyBucketIndicesY[7] += 1;
            }
            else
            {
                nearbyBucketIndicesY[2] -= 1;
                nearbyBucketIndicesY[3] -= 1;
                nearbyBucketIndicesY[6] -= 1;
                nearbyBucketIndicesY[7] -= 1;
            }

            if (( originIndex.z + 0.5f ) * cellSize <= position.z)
            {
                nearbyBucketIndicesZ[1] += 1;
                nearbyBucketIndicesZ[3] += 1;
                nearbyBucketIndicesZ[5] += 1;
                nearbyBucketIndicesZ[7] += 1;
            }
            else
            {
                nearbyBucketIndicesZ[1] -= 1;
                nearbyBucketIndicesZ[3] -= 1;
                nearbyBucketIndicesZ[5] -= 1;
                nearbyBucketIndicesZ[7] -= 1;
            }

            NativeArray<int> nearbyKeys = new NativeArray<int>(8, Allocator.Temp);
            for (int i = 0; i < 8; i++)
            {
                nearbyKeys[i] = SpatialHashing.Hash(new Vector3Int(nearbyBucketIndicesX[i], nearbyBucketIndicesY[i], nearbyBucketIndicesZ[i]), dimensions);
            }

            return nearbyKeys;
        }
    }
    
    [BurstCompile]
    private struct ComputeDensityPressure : IJobParallelFor
    {
        [ReadOnly] public NativeArray<int> neighbourTracker;
        [ReadOnly] public NativeArray<int> neighbourList; // Stores all neighbours of a particle aligned at 'particleIndex * maximumParticlesPerCell * 8'
        [ReadOnly] public NativeArray<Particle> particles;
        [NativeDisableParallelForRestriction] public NativeArray<float> densities;
        [NativeDisableParallelForRestriction] public NativeArray<float> pressures;
        
        [ReadOnly] public float mass;
        [ReadOnly] public float gasConstant;
        [ReadOnly] public float restDensity;
        [ReadOnly] public int maximumParticlesPerCell;
        [ReadOnly] public float radius2;
        [ReadOnly] public float radius3;
        
        public void Execute(int index)
        {
            // Doyub Kim 121, 122, 123
            // 5. Compute densities
            Vector3 origin = particles[index].Position;
            float sum = 0f;
            for (int j = 0; j < neighbourTracker[index]; j++)
            {
                int neighbourIndex = neighbourList[index * maximumParticlesPerCell * 8 + j];
                float distanceSquared = ( origin - particles[neighbourIndex].Position ).sqrMagnitude;
                sum += StdKernel(distanceSquared);
            }

            densities[index] = sum * mass + 0.000001f;

            // 6. Compute pressure based on density
            pressures[index] = gasConstant * ( densities[index] - restDensity ); // as described in Müller et al Equation 12
        }

        // Kernel by Müller et al.
        private float StdKernel(float distanceSquared)
        {
            // Doyub Kim
            float x = 1.0f - distanceSquared / radius2;
            return 315f / ( 64f * Mathf.PI * radius3 ) * x * x * x;
        }
    }

    [BurstCompile]
    private struct ComputeForces : IJobParallelFor
    {
        [ReadOnly] public NativeArray<int> neighbourTracker;
        [ReadOnly] public NativeArray<int> neighbourList; // Stores all neighbours of a particle aligned at 'particleIndex * maximumParticlesPerCell * 8'
        [ReadOnly] public NativeArray<Particle> particles;
        [NativeDisableParallelForRestriction] public NativeArray<Vector3> forces;
        [ReadOnly] public NativeArray<Vector3> velocities;
        [ReadOnly] public NativeArray<float> densities;
        [ReadOnly] public NativeArray<float> pressures;
        
        [ReadOnly] public float mass2;
        [ReadOnly] public int maximumParticlesPerCell;
        [ReadOnly] public float viscosityCoefficient;
        [ReadOnly] public Vector3 gravity;
        [ReadOnly] public float radius;
        [ReadOnly] public float radius4;
        [ReadOnly] public float radius5;
        
        public void Execute(int index)
        {
            forces[index] = Vector3.zero;
            var particleDensity2 = densities[index] * densities[index];
            for (int j = 0; j < neighbourTracker[index]; j++)
            {
                int neighbourIndex = neighbourList[index * maximumParticlesPerCell * 8 + j];
                float distance = ( particles[index].Position - particles[neighbourIndex].Position ).magnitude;
                if (distance > 0.0f)
                {
                    var direction = ( particles[index].Position - particles[neighbourIndex].Position ) / distance;
                    // 7. Compute pressure gradient force (Doyub Kim page 136)
                    forces[index] -= mass2 * ( pressures[index] / particleDensity2 + pressures[neighbourIndex] / ( densities[neighbourIndex] * densities[neighbourIndex] ) ) * SpikyKernelGradient(distance, direction);   // Kim
                    // 8. Compute the viscosity force
                    forces[index] += viscosityCoefficient * mass2 * ( velocities[neighbourIndex] - velocities[index] ) / densities[neighbourIndex] * SpikyKernelSecondDerivative(distance);    // Kim
                }
            }

            // Gravity
            forces[index] += gravity;
        }

        // Doyub Kim page 130
        private float SpikyKernelFirstDerivative(float distance)
        {
            float x = 1.0f - distance / radius;
            return -45.0f / ( Mathf.PI * radius4 ) * x * x;
        }

        // Doyub Kim page 130
        private float SpikyKernelSecondDerivative(float distance)
        {
            // Btw, it derives 'distance' not 'radius' (h)
            float x = 1.0f - distance / radius;
            return 90f / ( Mathf.PI * radius5 ) * x;
        }

        // // Doyub Kim page 130
        private Vector3 SpikyKernelGradient(float distance, Vector3 directionFromCenter)
        {
            return SpikyKernelFirstDerivative(distance) * directionFromCenter;
        }
    }

    [BurstCompile]
    private struct Integrate : IJobParallelFor
    {
        [NativeDisableParallelForRestriction] public NativeArray<Particle> particles;
        [ReadOnly] public NativeArray<Vector3> forces;
        [NativeDisableParallelForRestriction] public NativeArray<Vector3> velocities;

        [ReadOnly] public float mass;
        [ReadOnly] public float damping;
        [ReadOnly] public float dt;
        [ReadOnly] public int dimensions;
        
        public void Execute(int index)
        {
            var particle = particles[index];
            // forward Euler integration
            velocities[index] += dt * forces[index] / mass;
            particle.Position += dt * velocities[index];
            particles[index] = particle;
            
            particle = particles[index];
            var velocity = velocities[index];
            
            // enforce boundary conditions
            if (particles[index].Position.x - float.Epsilon < 0.0f)
            {
                velocity.x *= damping;
                particle.Position.x = float.Epsilon;
            }
            else if(particles[index].Position.x + float.Epsilon > dimensions - 1f) 
            {
                velocity.x *= damping;
                particle.Position.x = dimensions - 1 - float.Epsilon;
            }
            
            if (particles[index].Position.y - float.Epsilon < 0.0f)
            {
                velocity.y *= damping;
                particle.Position.y = float.Epsilon;
            }
            else if(particles[index].Position.y + float.Epsilon > dimensions - 1f) 
            {
                velocity.y *= damping;
                particle.Position.y = dimensions - 1 - float.Epsilon;
            }
            
            if (particles[index].Position.z - float.Epsilon < 0.0f)
            {
                velocity.z *= damping;
                particle.Position.z = float.Epsilon;
            }
            else if(particles[index].Position.z + float.Epsilon > dimensions - 1f) 
            {
                velocity.z *= damping;
                particle.Position.z = dimensions - 1 - float.Epsilon;
            }
            
            velocities[index] = velocity;
            particles[index] = particle;
        }
    }

    protected override void OnUpdate()
    {
        // Calculate hash of all particles and build neighboring list.
        // 1. Clear HashGrid
        _hashGrid.Clear();
        // 2. Recalculate hashes of each particle.
        RecalculateHashGrid recalculateHashGridJob = new RecalculateHashGrid
        {
            particles = _particles,
            hashGrid = _hashGrid.AsParallelWriter(),
            dimensions = SpatialHashing.Dimensions,
            cellSize = SpatialHashing.CellSize
        };
        JobHandle fillHashGridJobHandle = recalculateHashGridJob.Schedule(numberOfParticles, batchCount);
        fillHashGridJobHandle.Complete();
        // 3. For each particle go through all their 8 neighbouring cells.
        //    Check each particle in those neighbouring cells for interference radius r and store the interfering ones inside the particles neighbour list.
        BuildNeighbourList buildNeighbourListJob = new BuildNeighbourList
        {
            particles = _particles,
            hashGrid = _hashGrid,
            dimensions = SpatialHashing.Dimensions,
            cellSize = SpatialHashing.CellSize,
            radiusSquared = radius * radius,
            maximumParticlesPerCell = maximumParticlesPerCell,
            neighbourList = _neighbourList,
            neighbourTracker = _neighbourTracker
        };
        JobHandle buildNeighbourListJobHandle = buildNeighbourListJob.Schedule(numberOfParticles, batchCount);
        buildNeighbourListJobHandle.Complete();
        // 4. The Neighbouring-list should be n-particles big, each index containing a list of each particles neighbours in radius r.
        
        // 5. Compute density pressure
        ComputeDensityPressure computeDensityPressureJob = new ComputeDensityPressure
        {
            neighbourTracker = _neighbourTracker,
            neighbourList = _neighbourList,
            particles = _particles,
            densities = densities,
            pressures = pressures,
            mass = mass,
            gasConstant = gasConstant,
            restDensity = restDensity,
            maximumParticlesPerCell = maximumParticlesPerCell,
            radius2 = radius2,
            radius3 = radius3
        };
        JobHandle computeDensityPressureJobHandle = computeDensityPressureJob.Schedule(numberOfParticles, batchCount);
        computeDensityPressureJobHandle.Complete();

        ComputeForces computeForcesJob = new ComputeForces
        {
            neighbourTracker = _neighbourTracker,
            neighbourList = _neighbourList,
            particles = _particles,
            forces = forces,
            velocities = velocities,
            densities = densities,
            pressures = pressures,
            mass2 = mass2,
            maximumParticlesPerCell = maximumParticlesPerCell,
            viscosityCoefficient = viscosityCoefficient,
            gravity = gravity,
            radius = radius,
            radius4 = radius4,
            radius5 = radius5
        };
        JobHandle computeForcesJobHandle = computeForcesJob.Schedule(numberOfParticles, batchCount);
        computeForcesJobHandle.Complete();

        Integrate integrateJob = new Integrate
        {
            particles = _particles,
            forces = forces,
            velocities = velocities,
            mass = mass,
            damping = damping,
            dt = dt,
            dimensions = dimensions
        };
        JobHandle integrateJobHandle = integrateJob.Schedule(numberOfParticles, batchCount);
        integrateJobHandle.Complete();

        _particleColorPositionBuffer.SetData(_particles);
        material.SetFloat(SizeProperty, particleRenderSize);
        material.SetBuffer(ParticlesBufferProperty, _particleColorPositionBuffer);
        Graphics.DrawMeshInstancedIndirect(particleMesh, 0, material, new Bounds(Vector3.zero, new Vector3(100.0f, 100.0f, 100.0f)), _argsBuffer, castShadows: UnityEngine.Rendering.ShadowCastingMode.On);
        elapsedSimulationSteps++;
    }
}
