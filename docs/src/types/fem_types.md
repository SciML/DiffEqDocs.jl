# FEM Problems

Below are the definitions of the types which specify problems. Some general notes are:

* (t,x) vs (t,x,y): Mathematically one normally specifies equations in 2D as ``f(t,x,y)``.
  However, in this code we use `x` as a vector. Thus you can think of ``x``=`x[:,1]` and
  ``y``=`x[:,2]`. Thus input equations are of the form `f(x,t)` no matter the dimension.
  If time is not included in the problem (for example, a Poisson equation problem),
  then we use `f(x)`. An example is the equation ``u(x,y)= sin(2πx)cos(2πy)/(8π^2)``
  would be specified as `sol(x) = sin(2π.*x[:,1]).*cos(2π.*x[:,2])/(8π*π)`.
* Linearity: If the equation has a linear term, they are specified with functions
  `f(t,x)`. If it is nonlinear, it is specified with functions `f(t,x,u)`. The boundary
  conditions are always `(t,x)`
* Stochastic: By default the equation is deterministic. For each equation, one can
  specify a σ term which adds a stochastic ``σ(t,x,u)dW_t`` term to the equation
  (or with ``σ(t,x)dW_t`` if linear, must match `f`). ``dW_t`` corresponds to the type
  of noise which is chosen. By default this is space-time Gaussian white noise.

## Poisson Equation Problem

Wraps the data that defines a 2D linear Poisson equation problem:

```math
-Δu = f
```

with bounday conditions `gD` on the Dirichlet boundary and gN on the Neumann boundary.
Linearity is determined by whether the forcing function `f` is a function of one
variable `(x)` or two `(u,x)` (with `x=[:,1]` and `y=[:,2]`).

If the keyword `σ` is given, then this wraps the data that defines a 2D stochastic heat equation

```math
-Δu = f + σdW
```

### Constructors

`PoissonProblem(f,analytic,Du,mesh)`: Defines the Dirichlet problem with analytical solution `analytic`, solution gradient `Du = [u_x,u_y]`,
and forcing function `f`

`PoissonProblem(u0,f,mesh)`: Defines the problem with initial value `u0` (as a function) and f.
If your initial data is a vector, wrap it as `u0(x) = vector`.

Note: If all functions are of `(x)`, then the program assumes it's linear. Write
your functions using the math to program syntax translation: ``x`` `= x[:,1]` and ``y`` `= x[:,2]`.
Use `f=f(u,x)` and `σ=σ(u,x)` (if specified) for nonlinear problems
(with the boundary conditions still (x)). Systems of equations can be specified
with `u_i = u[:,i]` as the ith variable. See the example problems for more help.

### Keyword Arguments

* `gD` = Dirichlet boundary function

* `gN` = Neumann boundary function

* `σ` = The function which multiplies the noise ``dW``. By default `σ=0`.

* `noisetype` = A string which specifies the type of noise to be generated. By default
  `noisetype=:White` for Gaussian Spacetime White Noise.

* `numvars` = The number of variables in the Poisson system. Automatically calculated in many cases.

* `D` = Vector of diffusion coefficients. Defaults is `D=ones(1,numvars)`.

## Heat Equation Problem

Wraps the data that defines a 2D heat equation problem:

```math
u_t = Δu + f
```

with bounday conditions `gD` on the Dirichlet boundary and gN on the Neumann boundary.
Linearity is determined by whether the forcing function `f` is a function of two
variables `(t,x)` or three `(t,x,u)` (with `x=[:,1]` and `y=[:,2]`).

If the keyword `σ` is given, then this wraps the data that defines a 2D stochastic heat equation.

```math
u_t = Δu + f + σdW_t
```

### Constructors

* `HeatProblem(analytic,Du,f,mesh)`: Defines the Dirichlet problem with solution `analytic`,
  solution gradient `Du = [u_x,u_y]`, and the forcing function `f`.

* `HeatProblem(u0,f,mesh)`: Defines the problem with initial value `u0` (as a function) and `f`.
  If your initial data is a vector, wrap it as `u0(x) = vector`.

Note: If all functions are of `(t,x)`, then the program assumes it's linear. Write
your functions using the math to program syntax translation: ``x`` `= x[:,1]` and ``y`` `= x[:,2]`.
Use `f=f(t,x,u)` and `σ=σ(t,x,u)` (if specified) for nonlinear problems
(with the boundary conditions still (t,x)). Systems of equations can be specified
with `u_i = u[:,i]` as the ith variable. See the example problems for more help.

### Keyword Arguments

* `gD` = Dirichlet boundary function

* `gN` = Neumann boundary function

* `σ` = The function which multiplies the noise dW. By default `σ=0`.

* `noisetype` = A string which specifies the type of noise to be generated. By default
  `noisetype=:White` for Gaussian Spacetime White Noise.

* `numvars` = Number of variables in the system. Automatically calculated from u0 in most cases.

* `D` = Array which defines the diffusion coefficients. Default is `D=ones(1,numvars)`.

## Example Problems

Examples problems can be found in [DiffEqProblemLibrary.jl](https://github.com/JuliaDiffEq/DiffEqProblemLibrary.jl/blob/master/src/fem_premade_problems.jl).

To use a sample problem, you need to do:

```julia
# Pkg.add("DiffEqProblemLibrary")
using DiffEqProblemLibrary
```

### Poisson Equation

```@docs
DiffEqProblemLibrary.prob_poisson_birthdeathinteractingsystem
DiffEqProblemLibrary.prob_poisson_noisywave
DiffEqProblemLibrary.prob_poisson_birthdeathsystem
DiffEqProblemLibrary.prob_poisson_wave
DiffEqProblemLibrary.prob_poisson_birthdeath
```

### Heat Equation

```@docs
DiffEqProblemLibrary.prob_femheat_birthdeathsystem
DiffEqProblemLibrary.prob_femheat_birthdeathinteractingsystem
DiffEqProblemLibrary.prob_femheat_diffuse
DiffEqProblemLibrary.prob_femheat_stochasticbirthdeath
DiffEqProblemLibrary.prob_femheat_moving
DiffEqProblemLibrary.prob_femheat_pure
DiffEqProblemLibrary.prob_femheat_birthdeath
```
