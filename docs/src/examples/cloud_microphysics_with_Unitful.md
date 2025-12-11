# Atmospheric cloud microphysics in an adiabatic air-parcel model

This example constitutes a numerical implementation of a simple cloud physics model following :
**"An Elementary Parcel Model with Explicit Condensation and Supersaturation"**
by R.R. Rogers (1975).
The purpose of this examples is twofold: to provide an atmospheric-physics example for the package, and to demonstrate a robust way of using `DifferentialEquations.jl` with the [`Unitful.jl`](https://juliaphysics.github.io/Unitful.jl/stable/) dimensional analysis package.
This enables to programatically represent physical units in the code and to check dimensionality corectnes in all arithmetic opertions within physical equations, while incurring zero overhead within numerical solution.

The example code below reproduces all figures from the original paper by simulating the temporal evolution of relative humidity (i.e., saturation ratio), droplet radius, temperature, and liquid water content in an ascending air parcel.

* * *

## Physical Background

Clouds are suspensions of water droplets (and/or ice particles) in the air.
Formation of a cloud is a thermodynamic process.
To analyze it, we employ a simplified model based on the so-called air parcel framework, in which an adiabatically isolated "parcel" of air ascends  along a hydrostatic atmospheric pressure profile.
Adiabatic cooling of air due to its expansion causes increase of relative humidity, which upon reaching supersaturation, triggers diffusional growth of droplets and concurrent latent heat release associated with the vapor-liquid phase transition.

In this simplified model:

  - an **air parcel** rises at a constant velocity;
  - water vapor condenses on a fixed number of pre-existing monodisperse droplets (i.e., no aerosol interactions are modeled);
  - coalescence, sedimentation, radiative cooling and mixing are neglected.

Symbols used throughout the example below follow the notation from the Rogers 1975 paper.

* * *

## Model Setup

Initial conditions for the simulation include physical constants:

```@example rogers
import DifferentialEquations as DE
using Plots
using LaTeXStrings
using Unitful
using Test


consts_u = (
    R_prime = 287.053 * Unitful.u"J/(kg*K)",
    diffusion_consant = uconvert(Unitful.u"Pa*m^2/(K*s)", 8.28e2 *
                                                          Unitful.u"dyn/(K*s)"),
    g = 9.81 * Unitful.u"m/(s^2)",
    ε = 0.622 * Unitful.u"1",
    ρ_l = 1000.0 * Unitful.u"kg/(m^3)",
    cp = 1005 * Unitful.u"J/(kg*K)",
    L = 2.5e6 * Unitful.u"J/kg",
    convection_constant = uconvert(Unitful.u"W/(m*K)", 2.42e3 *Unitful.u"erg/(cm*s*K)"),
    # Equation A.1 and A.2 on page 204
    temp_273K = 273.0*Unitful.u"K",
    temp_393K = 393.0*Unitful.u"K",
    temp_120K = 120.0*Unitful.u"K",
    # Equation A.3 on page 204
    temp_eq_vapor = 5.44e3 * Unitful.u"K",
    eq_vapor_pressure = 2.75e12 * Unitful.u"mbar")

```

and following parameters:

  - Droplet concentration: ``n_\text{0}``
  - Updraft speed: ``U``
  - Initial droplet radius: ``r_\text{0}``
  - Initial temperature: ``T_\text{0}``
  - Initial pressure: ``p_\text{0}``
  - Initial saturation ratio: ``S_\text{0}`` 

They are defined by an initial state vector **u0** and a vector of parameters **params**.
In the paper, these values are used as representative of a rapidly developing cumulonimbus clouds in the lower troposphere.

The evolution of the state vector **u** is governed by the following equations:

  - Evolution of pressure in time (hydrostatic profile, eq. 5 on page 194):

```@example rogers
eqns = Dict{Symbol, Any}()
eqns[:p_time_rate] = (
    consts, eqns, params, p, T) -> -eqns.ρ(p, T) * consts.g * params.U

```

where ρ(p,T) is diagnosed from the ideal gas law:

```@example rogers
eqns[:ρ] = (consts, p, T) -> p / consts.R_prime / T
```

  - Evolution in time of droplet radius (diffusion of vapor in air, eq. 1 on page 193):

```@example rogers
eqns[:r_time_rate] = (consts, eqns, S, T, p, r) -> eqns.σ(S, T, p) / r
```

where σ is defined as (eq. 2 on page 193):

```@example rogers
eqns[:σ] = (
    consts, S, T, p) -> begin
    Fk, Fd = eqns_n.calculate_fk_fd(eqns_n, T, p)
    return (S - 1.0) / (Fk + Fd)
end
eqns[:calculate_fk_fd] = (consts,
    eqns,
    T,
    p) -> begin
    D = consts.diffusion_consant * T / p *
        (consts.temp_393K / (T + consts.temp_120K)) *
        (T / (consts.temp_273K))^(3/2)
    K = consts.convection_constant * (consts.temp_393K / (T + consts.temp_120K)) *
        (T / (consts.temp_273K))^(3/2)
    Fk = consts.L^2 * consts.ε * consts.ρ_l / (K * consts.R_prime * T^2)
    Fd = consts.R_prime * T * consts.ρ_l / (consts.ε * D * eqns.es(T))

    (Fk, Fd)
end
```

  - Evolution in time of liquid water mixing ratio (time derivative of eq. 4):

```@example rogers
eqns[:ksi_time_rate] = (consts, params, r,
    r_time_rate) -> 4.0 * π * consts.ρ_l * params.n_0 * r^2 * r_time_rate
```

where ``n_0`` is a constant parameter expressing the ratio of droplet number per mass of air.

  - Evolution in time of temperature (derived from eq. 6):

```@example rogers
eqns[:T_time_rate] = (consts,
    T,
    dp_dt,
    p,
    dksi_dt) -> T * consts.R_prime / consts.cp * dp_dt / p +
                consts.L / consts.cp * dksi_dt
```

  - Evolution of saturation ratio (eq. 10):

```@example rogers
eqns[:S_time_rate] = (consts, eqns, params, p, T,
    dksi_dt)->begin
    Q1, Q2 = eqns.calculate_q1_q2(eqns, T, p)
    return Q1 * params.U - eqns.ρ(p, T) * Q2 * dksi_dt
end
```

where coefficients Q1 and Q2 are calculated as:

```@example rogers
eqns[:calculate_q1_q2] = (consts,
    eqns,
    T,
    p) -> begin
    Q1 = consts.L * consts.g * consts.ε /
         (consts.R_prime * consts.cp * T^2) - consts.g / (consts.R_prime * T)
    Q2 = consts.R_prime * T / (consts.ε * eqns.es(T)) +
         consts.ε * consts.L^2 / (consts.cp * T * p)
    (Q1, Q2)
end
eqns[:es] = (
    consts, T) -> consts.eq_vapor_pressure * exp(-(consts.temp_eq_vapor) / T)
```

Addidionally we should define additional equations mentioned in the article such as proportionality factor G:

```@example rogers
eqns[:G] = (
    consts, eqns, T, p) -> begin
    Q1, Q2 = eqns.calculate_q1_q2(eqns, T, p)
    Fk, Fd = eqns.calculate_fk_fd(eqns, T, p)
    Q1 / Q2 * (Fk + Fd) / (4 * pi * consts.ρ_l * eqns.ρ(p, T))
end
```

which is used to define the limiting supersaturation (s=S-1, see eq. 16):

```@example rogers
eqns[:s_inf] = (consts, T, p, r, params) -> begin
    G_p_T = eqns_n.G(eqns_n, T, p)
    return G_p_T * params.U / r / params.n_0
end
```

As a first step in our calculation we will create a [NamedTuple](https://docs.julialang.org/en/v1/manual/types/#Named-Tuple-Types) from the dictionary 'eqns' declared above.

```@example rogers
eqns = (; eqns...)
```
### Testing dimensional corectness with Unitful.jl
We defined our constants using units from the article. To be sure our algorithm works well with units using Unitful.jl package. First, let's set first argument in our equations as constants with units.

```@example rogers
mapvalues(f, nt::NamedTuple) = NamedTuple{keys(nt)}(map(f, values(nt)))
eqns_u = mapvalues(f -> ((args...) -> f(consts_u, args...)), eqns)
```

Now we can test whether calculated density of air for some test values is within expected range.

```@example rogers
Test.@test 1 * Unitful.u"kg/(m^3)" <
           eqns_u.ρ(1000 * Unitful.u"hPa", 300 * Unitful.u"K") <
           1.5 * Unitful.u"kg/(m^3)"
```

As there is no error returned the test was successful.
Another test we can perform uses Figure 7 from the article depicting values of proportionality factor G in declared range.

```@example rogers
p_values = range(400 * Unitful.u"mbar", 1000 * Unitful.u"mbar", length = 50)
T_values = range(240 * Unitful.u"K", 300 * Unitful.u"K", length = 50)

Z = [uconvert(Unitful.u"s/g", eqns_u.G(eqns_u, T, p)*10^4)
     for p in p_values, T in T_values]
plt7 = Plots.contour(T_values, p_values, Z,
    xlabel = "T", ylabel = "p", yflip = true,
    title = "Function G(p,T) isolines")
plot(plt7)
```

Further calculations need to be performed without units as DifferentialEquations.jl solver expects vector with the initial state of the system without units.
We thus need to strip the constants from units and we can set it as the first argument in our equations (similarly as before).

```@example rogers
consts_n = map(Unitful.ustrip, consts_u)
eqns_n = mapvalues(f -> ((args...) -> f(consts_n, args...)), eqns)
```

Finally, it is possible to define a function that calculates the next step in time of state vector **u** using differential equations shown above.

```@example rogers
@enum vars p=1 T=2 S=3 r=4
function rhs!(du_dt, u, arguments, t)
    params = arguments.p
    f = arguments.f

    du_dt[Int(p)] = f.p_time_rate(f, params, u[Int(p)], u[Int(T)])

    du_dt[Int(r)] = f.r_time_rate(f, u[Int(S)], u[Int(T)], u[Int(p)], u[Int(r)])

    dksi_dt = f.ksi_time_rate(params, u[Int(r)], du_dt[Int(r)])

    du_dt[Int(T)] = f.T_time_rate(u[Int(T)], du_dt[Int(p)], u[Int(p)], dksi_dt)

    du_dt[Int(S)] = f.S_time_rate(f, params, u[Int(p)], u[Int(T)], dksi_dt)
end

```

To begin our simulation we have to define our initial conditions and parameters of the simulation.

```@example rogers
u0 = Vector{Quantity{Float64}}(undef, 4)
u0[Int(p)] = uconvert(Unitful.u"Pa", 800 * Unitful.u"mbar")
u0[Int(T)] = uconvert(Unitful.u"K", 7.0 * Unitful.u"°C")
u0[Int(S)] = 1.0 * Unitful.u"1"
u0[Int(r)] = uconvert(Unitful.u"m", 8 * Unitful.u"µm")

params = let
    U = 10Unitful.u"m/s"
    rho_0 = u0[Int(p)] / consts_u.R_prime / u0[Int(T)]
    n_0 = uconvert(Unitful.m^-3, 200Unitful.u"1/cm^3") / rho_0

    (U = U, rho_0 = rho_0, n_0 = n_0)
end
```

### Numerical solution using ODEProblem

Let's solve this system of equations using **ODEProblem** solver in a defined timespan.
To guarantee no integration-time overhead from using Unitful.jl, we strip the parameters and initial condition of units (the solver in fact demands it). 
We also pack parameters of the simulation and our equations into one NamedTuple.

```@example rogers
arguments = (
        p = NamedTuple{keys(params)}(Unitful.ustrip.(values(params))),
        f = eqns_n
    )
function simulate(;time_range=(0.0, 20.0), initial_saturation=1, parameters= params)
    initial_state = Unitful.ustrip.(u0)
    initial_state[Int(S)] = initial_saturation
    arguments = (
        p = NamedTuple{keys(parameters)}(Unitful.ustrip.(values(parameters))),
        f = eqns_n
    )

    return (DE.solve(DE.ODEProblem(rhs!, initial_state, time_range, arguments), saveat = 0.1), arguments.p)
end
solution,params = simulate()
```

### Figure 1 — Droplet Growth and Supersaturation

Now we can recreate Fig.1 from the article depicting profiles of pressure and radius.

```@example rogers
Plots.plot(solution.t, solution[Int(r), :] .* 1e6,
    ylabel = "Radius " * LaTeXStrings.L" [\mathrm{µm}]",
    xlabel = "Time t [s]",
    label = "r(t)")
Plots.plot!(twiny(), solution[Int(p), :] ./ 100, solution[Int(r), :] .* 1e6,
    xflip = true,
    label = false,
    xlabel = "Pressure [mbar]"
)
Plots.plot!(twinx(), solution.t, (solution[Int(S), :] .- 1) .* 100,
    ylabel = "Supersaturation [%]",
    color = :red,
    xticks = :none,
    label = "S(t)"
)
```

* * *

### Figure 2 — Liquid Water Content

Just as easily we can recreate Fig. 2 form the article.

```@example rogers
ksi = 4.0 * π * consts_n.ρ_l * params.n_0 * (solution[Int(r), :] .^ 3)

Plots.plot(solution.t, solution[Int(T), :],
    ylabel = "Temperature T [K]",
    xlabel = "Time t [s]",
    label = "T(t)",
    legend = :left
)

Plots.plot!(twiny(), solution[Int(p), :] ./ 100, solution[Int(T), :],
    xflip = true,
    label = false,
    xlabel = "Pressure [mbar]"
)
Plots.plot!(twinx(), solution.t, ksi,
    ylabel = "Liquid water mixing ratio "*LaTeXStrings.L" [\frac{\\mathrm{gm}}{\\mathrm{kg}}]",
    color = :red,
    xticks = :none,
    label = "x(t)",
    legend = :right
)
```

### Figure 3 and 4 - saturation ratio and radius initial conditions

To recreate Figre 3 and Figure 4, it is necessary to change the initial saturation ratio.

```@example rogers
plt3 = plot(title = "Supersaturation for different initial conditions",
    xlabel = "Time t [s]", ylabel = "S-1 [%]", framestyle = :box)

plt4 = plot(title = "Radius for different initial conditions", xlabel = "Time t [s]",
    ylabel = "Radius "*L" [\mathrm{µm}]", framestyle = :box)

for saturation in [1.02, 1.01, 1.0, 0.99, 0.98]
    solution2,_ = simulate(initial_saturation=saturation, time_range=(0, 30))
    Plots.plot!(plt3, solution2.t, (solution2[Int(S), :] .- 1)*100,
        label = "S(t=0) = $(round((saturation-1)*100,digits=2)) %"
    )
    Plots.plot!(plt4, solution2.t, solution2[Int(r), :] .* 1e6,
        label = "S(t=0) = $(round((saturation-1)*100, digits = 2)) %"
    )
end

plot(plt3)

```

```@example rogers
plot(plt4)
```


### Figure 5 - impact of the initial speed

For Figure 5, we need to change parameters of out simulation so we need to create new NamedTuple with our parameters each time.

```@example rogers
plt = Plots.plot(title = "Supersaturation for different updraft speed",
    xlabel = "Time t [s]",
    ylabel = "S-1 [%]"
)

for u in [1.0 * Unitful.u"m/s", 5.0 * Unitful.u"m/s", 15.0 * Unitful.u"m/s", 10.0*Unitful.u"m/s"]
    params = let
        U = u
        rho_0 = u0[Int(p)] / consts_u.R_prime / u0[Int(T)]
        n_0 = uconvert(Unitful.m^-3, 200Unitful.u"1/cm^3") / rho_0

        (U = U, rho_0 = rho_0, n_0 = n_0)
    end
    solution3,_ = simulate(parameters=params)
    Plots.plot!(plt, solution3.t, (solution3[Int(S), :] .- 1)*100,
        label = "U = $(u) "*LaTeXStrings.L" [\frac{\mathrm{m}}{\mathrm{s}}]"
    )
end

plot(plt)
```

### Figure 6 - impact of initial droplet concentration

In Figure 6, we need to change the parameter of initial concentration.

```@example rogers
plt = Plots.plot(title = "Supersaturation for various values \n of droplet concentra-
  tion",
    xlabel = "Time t [s]",
    ylabel = LaTeXStrings.L"\mathrm{ln}0(\frac{s_\infty-s}{s_0-s})",
    ylims = (-4, 0.5)
)

taus = Dict(50.0 => 5.45, 100.0 => 3.17, 500.0 => 0.73, 200.0 => 1.73)
colors = Dict(50.0 => "red", 100.0 => "green", 500.0 => "blue", 200.0 => "deeppink")

for n in [50.0, 100.0, 500.0, 200.0]
    params = let
        U = 10u"m/s"
        rho_0 = u0[Int(p)] / consts_u.R_prime / u0[Int(T)]
        n_0 = uconvert(Unitful.m^-3, n * Unitful.u"1/cm^3") / rho_0

        (U = U, rho_0 = rho_0, n_0 = n_0)
    end
    solution4,params = simulate(time_range=(0, 6), parameters=params)
    s_infty = eqns_n.s_inf.(solution4[Int(T), :], solution4[Int(p), :], solution4[Int(r), :], Ref(params))
    s=solution4[Int(S), :] .- 1
    line = (s_infty .- s) ./ (s_infty .- s[1])
    safe_line = max.(line, 0)
    line = log.(safe_line)
    Plots.plot!(plt, solution4.t, line,
        label = LaTeXStrings.L"n_0 = "*"$(n) "*LaTeXStrings.L" [\mathrm{cm}^{-3}]",
        color = colors[n]
    )
    exponential = -solution4.t ./ taus[n]
    Plots.plot!(plt, solution4.t, exponential;
        label = "",
        linestyle = :dash,
        color = colors[n]
    )
end
plot(plt)

```

* * *

## Reference

> R.R. Rogers (1975), *An Elementary Parcel Model with Explicit Condensation and Supersaturation*,
> Atmosphere, Vol. 13, No. 4, pp. 192–204.
> DOI: [10.1080/00046973.1975.9648397](https://doi.org/10.1080/00046973.1975.9648397)



