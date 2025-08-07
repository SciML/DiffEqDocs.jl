# Parcel Model with Explicit Condensation and Supersaturation

This example demonstrates a numerical implementation of a classic cloud physics model:  
**"An Elementary Parcel Model with Explicit Condensation and Supersaturation"**  
by R.R. Rogers (1975), using the Julia language and the [`DifferentialEquations.jl`](https://diffeq.sciml.ai/stable/) library.

The purpose of this example is to reproduce figures from the original paper by simulating the temporal evolution of supersaturation, droplet radius, temperature, and liquid water content in an ascending air parcel.

---

## Physical Background



Cloud formation is a complex thermodynamic process. To analyze it, we use simplified approaches, such as the  **parcel model**, in which a small "parcel" of air ascends adiabatically through the atmosphere.

In this simplified model:

- An **air parcel** rises at a constant velocity `U`.
- As pressure decreases, **supersaturation** is created.
- Water vapor condenses on pre-existing droplets (resulting in their growth), releasing latent heat.
- This process changes the **droplet radius**, **temperature**, and **liquid water content** over time.

Key simplifying assumptions:
- No new droplet activation (fixed number of droplets).
- No coalescence or sedimentation.
- No mixing with the environment (idealized adiabatic ascent).


---

##  Model Setup

Initial conditions for the simulation include physical constants: 

```@setup rogers
import DifferentialEquations as DE
using Plots
@enum vars p=1 T=2 S=3 r=4 

struct PhysicalConstants
    R::Float64    # J/(kg*K)
    g::Float64    # m/s^2
    eps::Float64
    ro_l::Float64 # kg/m^3
    cp::Float64   # J/(kg*K)
    L::Float64    # J/kg
    U::Float64    # m/s
end


const C = PhysicalConstants(
    287.053,     # R
    9.81,        # g
    0.622,       # eps
    1000.0,      # ro_l = 1 kg/L = 1000 kg/m^3
    1005.0,      # cp
    2.5e6,       # L = 2.5e3 J/g = 2.5e6 J/kg
    10.0         # U
)
```
and followin parameters:

- Initial droplet radius: $$r_0 $$
- Droplet concentration: $$n_0 $$
- Updraft speed: $$U $$
- Temperature: $$T_0 $$
- Pressure: $$p_0 $$
- Initial supersaturation: $$S_0 $$ (saturated)
- Gravitational acceleration $$g $$

```@example rogers

# Wektor początkowy u0
u0 = Vector{Float64}(undef, 4)
u0[Int(p)] = 800*100         # 800 mbar = 80000 Pa
u0[Int(T)] = 273.15 + 7.0    # 7°C = 280.15 K
u0[Int(S)] = 1.0             # bezwymiarowe
u0[Int(r)] = 8.0*10^(-6)          # 8 µm = 8e-6 m
param=zeros(2)
param[1]=10.0 # U
param[2]=u0[Int(p)]/C.R/u0[Int(T)]#rho_0
```
They are defined by an initial vector of a state function **u0** and a vector of parameters **param**.
These parameters are typical of rapidly developing cumulonimbus clouds in the lower troposphere.

The simulation integrates the following variables over time, which are described by following equations:
- Pressure 'p(t)' 
<div align="center">

$$
\frac{dp}{dt}  = -\rho_0 g U 
$$

</div>
where $$\rho = \frac{p}{R\cdot T}$$

- Droplet radius 'r(t)'
<div align="center">

$$
  \frac{dr}{dt} = \frac{\sigma(S, T, p)}{r}
$$

</div>
where $$\sigma = \frac{S-1}{F_D + F_K}$$


- Liquid water mixing ratio 'x(t)'
<div align="center">

$$
  \frac{dx}{dt} = 4.0 \cdot  \pi \cdot  \rho_l \cdot nu_0 \cdot r^2 \cdot \frac{dr}{dt}
$$

</div>
where $$nu_0 = 200\cdot10^(6) / \rho_0$$
- Temperature 'T(t)'
<div align="center">

$$
\frac{dT}{t}  = T \cdot \frac{R}{cp} \cdot \frac{\frac{dp}{dt} }{p} + \frac{L}{cp} \cdot  \frac{dx}{dt}
$$

</div>
- Supersaturation 'S(t)'
<div align="center">

$$
\frac{dS}{t}  = Q1 \cdot U - \rho \cdot Q2 \cdot \frac{dx}{dt}
$$

</div>
where $$Q1 = L \cdot g \cdot \frac{\epsilon}{R \cdot cp \cdot T^2 - \frac{g}{R \cdot T}}$$

Addidionally we should define all of funkctions and coefficients used in equations above.


```@example rogers
function es(T)
    e = 2.75 * 10^(12) * exp(-(5.44 * 10^3) / T)
    return e
end


function sigma(S, T, p)
    D = 2.6e-5#8.28*10^2 * T / p * (393.0 / (T + 120.0)) * (T / 273.0)^(3 / 2)
    println(D)
    K = 2.4e-2 #2.42*10^3 * (393.0 / (T + 120.0)) * (T / 273.0)^(3 / 2)
    println(K)
    Fk = C.L^2 * C.eps * C.ro_l / (K * C.R * T^2)
    Fd = C.R * T * C.ro_l / (C.eps * D * es(T))
    println(Fk," ", Fd)

    # println((S - 1.0) / (Fk + Fd))
    return (S - 1.0) / (Fk + Fd)
end


```

Finally, it is possible to define a function that calculates the next step in time of state function **u** using differenctial equations shown above. 

```@example rogers
function rhs!(du_dt, u, param, t)
    # eq 8
    ro = u[Int(p)] / C.R / u[Int(T)]
    # eq 5
    du_dt[Int(p)] = -ro * C.g * param[1]
    # # eq 1
    du_dt[Int(r)] = sigma(u[Int(S)], u[Int(T)], u[Int(p)]) / u[Int(r)]

    # rho0 i nu_0
    # rho0 = (80000.0) / (C.R * (7.0 + 273.15)) # 800 mbar -> 80000 Pa
    nu_0 = 200*10^(6) / param[2]  # 200 /cm^3 -> 2e6 /m^3

    # # eq 4
    dksi_dt = 4.0 * π * C.ro_l * nu_0 * u[Int(r)]^2 * du_dt[Int(r)]

    # # eq 6
    du_dt[Int(T)] = u[Int(T)] * C.R / C.cp * du_dt[Int(p)] / u[Int(p)] + C.L / C.cp * dksi_dt

    # # eq 10
    Q1 = C.L * C.g * C.eps / (C.R * C.cp * u[Int(T)]^2) - C.g / (C.R * u[Int(T)])
    # println(Q1)
    Q2 = C.R * u[Int(T)] / (C.eps * es(u[Int(T)])) + C.eps * C.L^2 / (C.cp * u[Int(T)] * u[Int(p)])
    # println(Q1 * param[1], - ro * Q2 * dksi_dt)
    du_dt[Int(S)] = Q1 * param[1] - ro * Q2 * dksi_dt
end


```
Let's solve this system of equations using **ODEProblem** solver in defined timespan. 

```@example rogers
tspan = (0.0, 20.0)
prob = DE.ODEProblem(rhs!, u0, tspan,param)
sol = DE.solve(prob)

```


---

## Results

Using the solution from above  it is possible to recreate figures from the original paper.

### Figure 1 — Droplet Growth and Supersaturation
```@example rogers
Plots.plot(sol.t, sol[Int(r),:] .* 1e6, 
    ylabel="Radius [um]", xlabel="Time [s]", label = "r(t)"
)
Plots.plot!(twinx(),sol.t, (sol[Int(S),:] .- 1).*100,
    ylabel="Supersaturation [%]",
    color=:red,
    xticks=:none,
    label="S(t)"
)
savefig( "rogers_fig1.svg")
```



---

### Figure 2 — Temperature and Liquid Water Content
```@example rogers
nu_0 = 200*10^(6) / param[2]  # 200 /cm^3 -> 2e6 /m^3
ksi  =  4.0 * π * C.ro_l * nu_0 *(sol[Int(r),:].^3)
Plots.plot(sol.t, sol[Int(T),:],
    ylabel="Temperatura [K]", xlabel="Czas [s]", label = "T(t)")
Plots.plot!(twinx(),sol.t, ksi,
    ylabel="Współczynnik mieszania [gm/kg]",
    color=:red,
    xticks=:none,
    label="x(t)"
)
savefig( "rogers_fig2.svg")

```


---

## Reference

> R.R. Rogers (1975), *An Elementary Parcel Model with Explicit Condensation and Supersaturation*,  
> Atmosphere, Vol. 13, No. 4, pp. 192–204.  
> DOI: [10.1080/00046973.1975.9648397](https://doi.org/10.1080/00046973.1975.9648397)

