# [Noise Processes](@id noise_process)

Noise processes are essential in continuous stochastic modeling. The `NoiseProcess`
types are distributionally-exact, meaning they are not solutions of
stochastic differential equations and instead are directly generated according
to their analytical distributions. These processes are used as the noise term
in the SDE and RODE solvers. Additionally, the noise processes themselves can
be simulated and solved using the DiffEq common interface (including the Monte
Carlo interface).

For more details, see [DiffEqNoiseProcess.jl](noise.sciml.ai/dev)


