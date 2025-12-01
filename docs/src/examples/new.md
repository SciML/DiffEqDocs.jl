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
using LaTeXStrings
using Unitful
using Test
@enum vars p=1 T=2 S=3 r=4 

constants_u = (
    R_prime = 287.053 * Unitful.u"J/(kg*K)",
    diffusion_consant =  uconvert(Unitful.u"Pa*m^2/(K*s)", 8.28*10^2 * Unitful.u"dyn/(K*s)"),
    g = 9.81 * Unitful.u"m/(s^2)",
    ε = 0.622 * Unitful.u"1",
    ρ_l = 1000.0 * Unitful.u"kg/(m^3)",
    cp = 1005 * Unitful.u"J/(kg*K)",
    L = 2.5e6 * Unitful.u"J/kg",
    convection_constant = uconvert(Unitful.u"W/(m*K)", 2.42*10^3 * Unitful.u"erg/(cm*s*K)"),
    temp_273K = 273.0*Unitful.u"K",
    temp_393K = 393.0*Unitful.u"K",
    temp_120K = 120.0*Unitful.u"K",
    temp_eq_vapor = 5.44 * 10^3 *Unitful.u"K", 
    eq_vapor_pressure = 2.75 * 10^(12)*Unitful.u"mbar"

)

```
and following parameters:

- Initial droplet radius: ``r_\text{0} ``
- Droplet concentration: ``n_\text{0} ``
- Updraft speed: ``U ``
- Temperature: ``T_\text{0} ``
- Pressure: ``p_\text{0} ``
- Initial supersaturation: ``S_\text{0} `` (saturated)


They are defined by an initial vector of a state function **u0** and a vector of parameters **param**.
These parameters are typical of rapidly developing cumulonimbus clouds in the lower troposphere.

The simulation integrates the following variables over time, which are described by following equations:
- Pressure 'p(t)' 


```@example rogers
formulas = Dict{Symbol, Any}()
formulas[:p_time_rate] = (constants, formulas, param, p, T) -> -formulas.ρ(p,T) * constants.g * param.U

```
where ρ is a function representing density:
```@example rogers
formulas[:ρ] = (constants, p, T) -> p / constants.R_prime / T
```






- Droplet radius 'r(t)'

```@example rogers
formulas[:r_time_rate] = (constants,formulas, S,T,p,r) -> formulas.σ(S, T, p) / r
```
where σ is a function
```@example rogers
formulas[:σ] = (constants,S, T, p)-> begin 
        Fk, Fd = formulas_n.calculate_fk_fd(formulas_n,T,p)
        return (S - 1.0) / (Fk + Fd)
    end
formulas[:calculate_fk_fd]= (constants,formulas, T, p) -> begin
        D = constants.diffusion_consant * T / p * (constants.temp_393K / (T + constants.temp_120K)) * (T / (constants.temp_273K))^(3/2)
        K = constants.convection_constant * (constants.temp_393K / (T + constants.temp_120K)) * (T / (constants.temp_273K))^(3/2)
        Fk = constants.L^2 * constants.ε * constants.ρ_l / (K * constants.R_prime * T^2)
        Fd = constants.R_prime * T * constants.ρ_l / (constants.ε * D * formulas.es(T))

        (Fk, Fd)
    end
```




- Liquid water mixing ratio 'x(t)'


```@example rogers
formulas[:ksi_time_rate] = (constants, param, r, r_time_rate) -> 4.0 * π * constants.ρ_l * param.n_0 * r^2 * r_time_rate
```

where ``nu_0 `` is a parameter set for the simulation
- Temperature 'T(t)'


```@example rogers
formulas[:T_time_rate] = (constants, T, dp_dt, p, dksi_dt) -> T * constants.R_prime / constants.cp * dp_dt / p + constants.L / constants.cp * dksi_dt
```


- Supersaturation 'S(t)'


```@example rogers
formulas[:S_time_rate] = (constants,formulas,param, p, T, dksi_dt)->begin
        Q1,Q2 = formulas.calculate_q1_q2(formulas,T,p)
        return Q1 * param.U - formulas.ρ(p,T) * Q2 * dksi_dt
    end
```
where coefficients Q1 and Q2 are calculated as:
```@example rogers
formulas[:calculate_q1_q2] = (constants,formulas,T, p) -> begin
      
      Q1 = constants.L * constants.g * constants.ε / (constants.R_prime * constants.cp * T^2) - constants.g / (constants.R_prime * T)
      Q2 = constants.R_prime * T/ (constants.ε * formulas.es(T)) + constants.ε * constants.L^2 / (constants.cp * T * p)
      (Q1,Q2)
    end
formulas[:es] = (constants,T) -> constants.eq_vapor_pressure * exp(-(constants.temp_eq_vapor) / T)
```


Addidionally we should define additional formulas mentioned in the article such as proportionality factor
```@example rogers
formulas[:G] = (constants,formulas, T, p) -> begin
        Q1,Q2 = formulas.calculate_q1_q2(formulas,T,p)
        Fk, Fd  = formulas.calculate_fk_fd(formulas,T,p)
        Q1/Q2*(Fk+Fd)/(4*pi*constants.ρ_l*formulas.ρ(p,T))
    end
```

and limiting supersaturation
```@example rogers
formulas[:s_inf] = (constants,T,p, r,param) -> begin
        G_p_T = formulas_n.G(formulas_n,T,p)
        return G_p_T*param.U/r/param.n_0
    end
```
As a first step in our calculation we will create a NamedTuple from the dictionary 'formulas' declared above.
```@example rogers
formulas = (; formulas...)
```
We defined our constants using units from the article. To be sure our algorith works well with units using Unitful.jl package. First, let's set first argument in our formulas as constants with units. 
```@example rogers
mapvalues(f, nt::NamedTuple) = NamedTuple{keys(nt)}(map(f, values(nt)))
formulas_u = mapvalues(f -> ((args...) -> f(constants_u, args...)), formulas)   
```
Now we can test whether calulated density of air for some test values is within expected range.
```@example rogers
@Test.test 1 * Unitful.u"kg/(m^3)" < formulas_u.ρ(1000 * Unitful.u"hPa", 300 * Unitful.u"K") < 1.5 * Unitful.u"kg/(m^3)"
```
As there is no error returned the test was succesful. 
Another test we can perform uses Figure 7 from the article depicting values of proportionality factor G in declared range.
```@example rogers
p_values = range(400* Unitful.u"mbar", 1000 * Unitful.u"mbar", length=50)
T_values = range(240 * Unitful.u"K", 300 * Unitful.u"K", length=50)

Z = [uconvert(Unitful.u"s/g",formulas_u.G(formulas_u,T,p)*10^4) for p in p_values, T in T_values] 
plt7 =Plots.contour( T_values,p_values, Z, 
        xlabel="T", ylabel="p", yflip=true,
        title="Function G(p,T) isolines")
plot(plt7)
```

Further calculations need to be performed without units as DifferentialEquations.jl solver expects vector with the begigning state of the system without units. We need to strip constants form units and we can set it as the first argument in our formulas (similarly as before).
```@example rogers
constants_n = map(Unitful.ustrip, constants_u)
formulas_n = mapvalues(f -> ((args...) -> f(constants_n, args...)), formulas)   
```


Finally, it is possible to define a function that calculates the next step in time of state function **u** using differenctial equations shown above. 

```@example rogers
    function rhs!(du_dt, u, parameters, t)

    param = parameters.p
    f = parameters.f

    du_dt[Int(p)] = f.p_time_rate(f, param,u[Int(p)],u[Int(T)])

    du_dt[Int(r)] = f.r_time_rate(f, u[Int(S)],u[Int(T)],u[Int(p)],u[Int(r)])
    
    dksi_dt =  f.ksi_time_rate( param, u[Int(r)], du_dt[Int(r)])

    du_dt[Int(T)] = f.T_time_rate( u[Int(T)], du_dt[Int(p)], u[Int(p)], dksi_dt)


    du_dt[Int(S)] = f.S_time_rate(f,param, u[Int(p)], u[Int(T)], dksi_dt)

end


```
To begin our simulation we have to define our begining conditions and parameters of the simulation. 
```@example rogers
u0 = Vector{Quantity{Float64}}(undef, 4)
u0[Int(p)] = uconvert(Unitful.u"Pa", 800 * Unitful.u"mbar")       
u0[Int(T)] = uconvert(Unitful.u"K", 7.0 * Unitful.u"°C")  
u0[Int(S)] = 1.0 * Unitful.u"1"             
u0[Int(r)] = uconvert(Unitful.u"m", 8 * Unitful.u"µm")
       # 
param = let
    U = 10Unitful.u"m/s"
    rho_0 = u0[Int(p)] / constants_u.R_prime / u0[Int(T)]
    n_0 = uconvert(Unitful.m^-3, 200Unitful.u"1/cm^3") / rho_0

    (U = U, rho_0 = rho_0, n_0 = n_0)
end
```
Next we have to strip them of units (as solver demands). We also pack parameters of the simulation and our formulas into one NamedTuple. 
```@example rogers
u0 = Unitful.ustrip.(u0)
param = NamedTuple{keys(param)}(Unitful.ustrip.(values(param)))

parameters = (
    p = param,
    f = formulas_n,
)

```
Let's solve this system of equations using **ODEProblem** solver in defined timespan. 

```@example rogers
solution = DE.solve(DE.ODEProblem(rhs!, u0, (0.0, 20.0),parameters), saveat = 0.1)


```
### Figure 1 — Droplet Growth and Supersaturation
Now we can recreate Fig.1 from the article depiction pressure and radius when our parcel of air is rising.
```@example rogers
Plots.plot(solution.t, solution[Int(r),:] .* 1e6, 
        ylabel="Radius " * LaTeXStrings.L" [\mathrm{µm}]", 
        xlabel="Time [s]", 
        label = "r(t)")
Plots.plot!(twiny(), solution[Int(p),:]./100, solution[Int(r),:] .* 1e6, 
        xflip = true, 
        label = false,
        xlabel = "Pressure [mbar]" 
    )
Plots.plot!(twinx(),solution.t, (solution[Int(S),:] .- 1).*100,
        ylabel="Supersaturation [%]",
        color=:red,
        xticks=:none,
        label="S(t)"
    )
```



---

### Figure 2 — Liquid Water Content
Just as easily we can recreate Fig. 2 form the article. 
```@example rogers
ksi  =  4.0 * π * constants_n.ρ_l * parameters.p.n_0 *(solution[Int(r),:].^3)

Plots.plot(solution.t, solution[Int(T),:],
        ylabel="Temperature [K]", 
        xlabel="Time [s]", 
        label = "T(t)",
        legend=:left
    )

Plots.plot!(twiny(),solution[Int(p),:]./100, solution[Int(T),:], 
        xflip = true, 
        label = false,
        xlabel = "Pressure [mbar]" 
    )
Plots.plot!(twinx(),solution.t, ksi,
        ylabel="Liquid water mixing ratio "*LaTeXStrings.L" [\frac{\\mathrm{gm}}{\\mathrm{kg}}]",
        color=:red,
        xticks=:none,
        label="x(t)",
        legend=:right
    )
```

### Figure 3 and 4 - supersaturation and radius initial conditions
To recreate Figre 3 and Figure 4 it is necessary to change initial condition of supersaturation for the simulation.

```@example rogers
plt3 = plot(title="Supersaturation for different initial conditions", xlabel="t", ylabel="S-1 [%]", framestyle = :box)

plt4 = plot(title="Radius for different initial conditions", xlabel="t", ylabel="Radius "*L" [\mathrm{µm}]", framestyle = :box)

for s in [1.02,1.01,0.99, 0.98,1.0]
    u0[Int(S)] = s
    solution2 = DE.solve(DE.ODEProblem(rhs!, u0, (0.0, 30.0),parameters), saveat = 0.1)
    Plots.plot!(plt3, solution2.t, (solution2[Int(S),:] .- 1)*100, 
    label="S(t=0) = $(round((s-1)*100,digits=2)) %"
    )
    Plots.plot!(plt4, solution2.t, solution2[Int(r),:] .* 1e6, 
    label="S(t=0) = $(round((s-1)*100, digits = 2)) %"
    )

end

plot(plt3)   

```

```@example rogers
plot(plt4)
u0[Int(S)]=1.0
```
At the end we return to the previous value of initial supersaturation.

### Figure 5 - impact of the initial speed
For Figure 5 we need to change parameters of out simulation so we need to create new NamedTuple with our parameters each time. 

```@example rogers
plt = Plots.plot(title="Supersaturation for different updraft speed", 
    xlabel="Time T[s]", 
    ylabel="S-1 [%]"
)


for u in [1.0 * Unitful.u"m/s",5.0* Unitful.u"m/s",15.0* Unitful.u"m/s",10.0*Unitful.u"m/s"]
    param = let
        U = u 
        rho_0 = u0[Int(p)] / constants_u.R_prime / u0[Int(T)]
        n_0 = uconvert(Unitful.m^-3, 200Unitful.u"1/cm^3") / rho_0

        (U = U, rho_0 = rho_0, n_0 = n_0)
    end
    param = NamedTuple{keys(param)}(Unitful.ustrip.(values(param)))

    parameters = (
        p = param,
        f = formulas_n,
    )
    solution3 = DE.solve(DE.ODEProblem(rhs!, u0, (0.0, 20.0),parameters), saveat = 0.1)
    Plots.plot!(plt, solution3.t, (solution3[Int(S),:] .- 1)*100, 
    label="U = $(u) "*LaTeXStrings.L" [\frac{\mathrm{m}}{\mathrm{s}}]"
    )


end

plot(plt)
```

### Figure 6 - impact of initial droplet concentration
In Figure 6 we need to change the parameter of initial concentration. 
```@example rogers
plt = Plots.plot(title="Supersaturation for various values \n of droplet concentra-
tion", 
    xlabel="t", 
    ylabel=LaTeXStrings.L"\mathrm{ln}0(\frac{s_\infty-s}{s_0-s})",
    ylims = (-4,0.5)
)

taus = Dict(50.0 => 5.45, 100.0 => 3.17,500.0 => 0.73,200.0 => 1.73)
colors = Dict(50.0 => "red", 100.0 => "green",500.0 => "blue",200.0 => "deeppink")

for n in [50.0,100.0,500.0,200.0]
    param = let
        U = 10u"m/s"
        rho_0 = u0[Int(p)] / constants_u.R_prime / u0[Int(T)]
        n_0 = uconvert(Unitful.m^-3, n * Unitful.u"1/cm^3") / rho_0

        (U = U, rho_0 = rho_0, n_0 = n_0)
    end
    param = NamedTuple{keys(param)}(Unitful.ustrip.(values(param)))

    parameters = (
        p = param,
        f = formulas_n,
    )
    solution4 = DE.solve(DE.ODEProblem(rhs!, u0, (0.0,6.0) ,parameters), saveat = 0.1)
    s_infty = formulas_n.s_inf.(solution4[Int(T),:],solution4[Int(p),:],solution4[Int(r),:], Ref(param))
    s=solution4[Int(S),:] .- 1
    line = (s_infty.-s)./(s_infty.-s[1])
    safe_line = max.(line, 0)
    line = log.(safe_line)
        Plots.plot!(plt, solution4.t, line, 
    label=LaTeXStrings.L"n_0 = "*"$(n) "*LaTeXStrings.L" [\mathrm{cm}^{-3}]",
    color = colors[n]
    )
    exponential = -solution4.t./taus[n]
    Plots.plot!(plt, solution4.t, exponential;
    label = "",         
    linestyle = :dash,   
    color = colors[n]        
    )



end
plot(plt)



```


---

## Reference

> R.R. Rogers (1975), *An Elementary Parcel Model with Explicit Condensation and Supersaturation*,  
> Atmosphere, Vol. 13, No. 4, pp. 192–204.  
> DOI: [10.1080/00046973.1975.9648397](https://doi.org/10.1080/00046973.1975.9648397)




