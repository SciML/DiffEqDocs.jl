# Plot Functions

## Standard Plots

Plotting functionality is provided by recipes to Plots.jl. To
use plot solutions, simply call the `plot(type)` after importing Plots.jl
and the plotter will generate appropriate plots.

```julia
using Plots
plot(sol) # Plots the solution
```

Many of the types defined in the DiffEq universe, such as
`ODESolution`, `ConvergenceSimulation` `WorkPrecision`, etc. have plot recipes
to handle the default plotting behavior. Plots can be customized using
[all of the keyword arguments provided by Plots.jl](https://juliaplots.github.io/supported/).
For example, we can change the plotting backend to the GR package and put a title
on the plot by doing:

```julia
gr()
plot(sol,title="I Love DiffEqs!")
```

## Animations

Using the iterator interface over the solutions, animations can also be generated
via the `animate(sol)` command. One can choose the `filename` to save to,
frames per second `fps`, and the density of steps to show `every` via keyword arguments.
The rest of the arguments will be directly passed to the plot recipe to be handled
as normal. For example, we can animate our solution with a larger line-width which
saves every 4th frame via:

```julia
animate(sol,lw=3,every=4)
```

Please see [Plots.jl's documentation](https://juliaplots.github.io/) for more information
on the available attributes.
