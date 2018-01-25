# Physical Models

The physical modeling functionality is provided by DiffEqPhysics.jl and helps
the user build and solve the differential equation based physical models.

## Hamiltonian Problems

ODEs defined by Hamiltonians is described in the
[Dynamical ODEs section](../../types/dynamical_types.html).

## N-Body Problems

```julia
nprob = NBodyProblem(f, mass, vel, pos, tspan)
```

where `f` is the potential function, `mass` is the mass matrix, `pos` and `vel`
are `ArrayPartition`s for the intial positions and velocity, and `tspan` is the
timespan to solve on.

### Example

In this example we will model the outer solar system planets.

```julia
G = 2.95912208286e-4
M = [1.00000597682, 0.000954786104043, 0.000285583733151, 0.0000437273164546, 0.0000517759138449, 1/1.3e8]
invM = inv.(M)
planets = ["Sun", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"]

pos_x = [0.0,-3.5023653,9.0755314,8.3101420,11.4707666,-15.5387357]
pos_y = [0.0,-3.8169847,-3.0458353,-16.2901086,-25.7294829,-25.2225594]
pos_z = [0.0,-1.5507963,-1.6483708,-7.2521278,-10.8169456,-3.1902382]
pos = ArrayPartition(pos_x,pos_y,pos_z)

vel_x = [0.0,0.00565429,0.00168318,0.00354178,0.00288930,0.00276725]
vel_y = [0.0,-0.00412490,0.00483525,0.00137102,0.00114527,-0.00170702]
vel_z = [0.0,-0.00190589,0.00192462,0.00055029,0.00039677,-0.00136504]
vel = ArrayPartition(vel_x,vel_y,vel_z)

tspan = (0.,200_000)

const ∑ = sum
const N = 6
potential(p, t, x, y, z, M) = -G*∑(i->∑(j->(M[i]*M[j])/sqrt((x[i]-x[j])^2 + (y[i]-y[j])^2 + (z[i]-z[j])^2), 1:i-1), 2:N)
nprob = NBodyProblem(potential, M, vel, pos, tspan)
sol = solve(nprob,Yoshida6(), dt=100)
```
