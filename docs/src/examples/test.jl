using DifferentialEquations
using Unitful
using Test
using LaTeXStrings
using Plots
@enum vars p=1 T=2 S=3 r=4 
@enum par U=1 rho_0=2 n_0=3



#functions using constants
functions = (
    rho = (constants, pressure, temperature) -> pressure / constants.specific_gas_constant / temperature,
    calculate_fk_fd = (constants,functions, temperature, pressure) -> begin
        D = constants.diffusion_consant * temperature / pressure * (constants.temp_393K / (temperature + constants.temp_120K)) * (temperature / (constants.temp_273K))^(3/2)
        K = constants.convection_constant * (constants.temp_393K / (temperature + constants.temp_120K)) * (temperature / (constants.temp_273K))^(3/2)
        Fk = constants.L^2 * constants.epsilon * constants.rho_l / (K * constants.specific_gas_constant * temperature^2)
        Fd = constants.specific_gas_constant * temperature * constants.rho_l / (constants.epsilon * D * functions.es(temperature))

        (Fk, Fd)
    end,
    calculate_q1_q2 = (constants,functions,temperature, pressure) -> begin
      
      Q1 = constants.L * constants.g * constants.epsilon / (constants.specific_gas_constant * constants.cp * temperature^2) - constants.g / (constants.specific_gas_constant * temperature)
      Q2 = constants.specific_gas_constant * temperature/ (constants.epsilon * functions.es(temperature)) + constants.epsilon * constants.L^2 / (constants.cp * temperature * pressure)
      (Q1,Q2)
    end,
    G = (constants,functions, temperature, pressure) -> begin
        Q1,Q2 = functions.calculate_q1_q2(functions,temperature,pressure)
        Fk, Fd  = functions.calculate_fk_fd(functions,temperature,pressure)
        Q1/Q2*(Fk+Fd)/(4*pi*constants.rho_l*functions.rho(pressure,temperature))
    end, 
    es = (constants,temperature) -> constants.eq_vapor_pressure * exp(-(constants.temp_eq_vapor) / temperature)
    )

#physical constants
constants_u = (
    specific_gas_constant = 287.053 * u"J/(kg*K)",
    diffusion_consant =  uconvert(u"Pa*m^2/(K*s)", 8.28*10^2 * u"dyn/(K*s)"),
    g = 9.81 * u"m/(s^2)",
    epsilon = 0.622 * u"1",
    rho_l = 1000.0 * u"kg/(m^3)",
    cp = 1005 * u"J/(kg*K)",
    L = 2.5e6 * u"J/kg",
    convection_constant = uconvert(u"W/(m*K)", 2.42*10^3 * u"erg/(cm*s*K)"),
    temp_273K = 273.0*u"K",
    temp_393K = 393.0*u"K",
    temp_120K = 120.0*u"K",
    temp_eq_vapor = 5.44 * 10^3 *u"K", 
    eq_vapor_pressure = 2.75 * 10^(12)*u"mbar"

)






#testing wether units are correct
begin 
    mapvalues(f, nt::NamedTuple) = NamedTuple{keys(nt)}(map(f, values(nt)))
    functions_u = mapvalues(f -> ((args...) -> f(constants_u, args...)), functions)    
    @test 1 * u"kg/(m^3)" < functions_u.rho(1000 * u"hPa", 300 * u"K") < 1.5 * u"kg/(m^3)"
    @show Unitful.dimension(constants_u.g)

    p_min = 400* u"mbar"
    p_max = 1000 * u"mbar"
    T_min = 240 * u"K"
    T_max = 300 * u"K"
    p_values = range(p_min, p_max, length=50)
    T_values = range(T_min, T_max, length=50)

    #tu jest blad w wartosciach 
    Z = [uconvert(u"s/g",functions_u.G(functions_u,T,p)) for p in p_values, T in T_values] 
    plt =Plots.contour( T_values,p_values, Z, 
            xlabel="T", ylabel="p", yflip=true,
            title="Izolinie funkcji G(p,T)")
    savefig(plt,"rogers_fig7.svg")
end


#calculations without units
begin
    constants_n = map(ustrip, constants_u)
    functions_n = mapvalues(f -> ((args...) -> f(constants_n, args...)), functions)   

    #functions not using constants
    function sigma(S, T, p)
        Fk, Fd = functions_n.calculate_fk_fd(functions_n,T,p)
        return (S - 1.0) / (Fk + Fd)
    end

    function s_inf(T,p, r,param)
        G_p_T = functions_n.G(functions_n,T,p)
        return G_p_T*param[Int(U)]/r/param[Int(n_0)]
    end

    function rhs!(du_dt, u, param, t)
   
        # eq 5
        du_dt[Int(p)] = -functions_n.rho(u[Int(p)],u[Int(T)]) * constants_n.g * param[1]
        # eq 1
        du_dt[Int(r)] = sigma(u[Int(S)], u[Int(T)], u[Int(p)]) / u[Int(r)]
        
        # eq 4
        dksi_dt = 4.0 * π * constants_n.rho_l * param[Int(n_0)] * u[Int(r)]^2 * du_dt[Int(r)]

        # eq 6
        du_dt[Int(T)] = u[Int(T)] * constants_n.specific_gas_constant / constants_n.cp * du_dt[Int(p)] / u[Int(p)] + constants_n.L / constants_n.cp * dksi_dt

        # eq 10
        Q1,Q2 = functions_n.calculate_q1_q2(functions_n,u[Int(T)],u[Int(p)])

        du_dt[Int(S)] = Q1 * param[1] - functions_n.rho(u[Int(p)],u[Int(T)]) * Q2 * dksi_dt
    end

    # Wektor początkowy u0
    u0 = Vector{Quantity{Float64}}(undef, 4)
    u0[Int(p)] = uconvert(u"Pa", 800 * u"mbar")       
    u0[Int(T)] = uconvert(u"K", 7.0 * u"°C")  
    u0[Int(S)] = 1.0 * u"1"             # bezwymiarowe
    u0[Int(r)] = uconvert(u"m", 8 * u"µm")          # 8 µm = 8e-6 m

    #parameters
    param=Vector{Quantity}(undef, 3)
    param[Int(U)]=10.0 * u"m/s"
    param[Int(rho_0)]=u0[Int(p)]/constants_u.specific_gas_constant/u0[Int(T)]
    param[Int(n_0)] = uconvert(Unitful.m^-3, 200.0u"1/cm^3") / param[Int(rho_0)]


    u0 = ustrip.(u0)
    param = ustrip.(param)
    # Zakres czasu 
    tspan = (0.0, 20.0)

    # Definicja i rozwiązanie problemu ODE
    prob = ODEProblem(rhs!, u0, tspan,param)
    sol = solve(prob, saveat = 0.1)

    # Wykres r i S
    plot(sol.t, sol[Int(r),:] .* 1e6, 
        ylabel="Radius " * L" [\mathrm{µm}]", 
        xlabel="Time [s]", 
        label = "r(t)")
    plot!(twiny(),sol[Int(p),:]./100, sol[Int(r),:] .* 1e6, 
        xflip = true, 
        label = false,
        xlabel = "Pressure [mbar]" 
    )
    plot!(twinx(),sol.t, (sol[Int(S),:] .- 1).*100,
        ylabel="Supersaturation [%]",
        color=:red,
        xticks=:none,
        label="S(t)"
    )
    savefig( "rogers_fig1.svg")


    #Wykres ksi i T
    ksi  =  4.0 * π * constants_n.rho_l * param[Int(n_0)] *(sol[Int(r),:].^3)
    plot(sol.t, sol[Int(T),:],
        ylabel="Temperature [K]", 
        xlabel="Time [s]", 
        label = "T(t)",
        legend=:left
    )

    plot!(twiny(),sol[Int(p),:]./100, sol[Int(T),:], 
        xflip = true, 
        label = false,
        xlabel = "Pressure [mbar]" 
    )
    plot!(twinx(),sol.t, ksi,
        ylabel="Liquid water mixing ratio "*L" [\frac{\mathrm{gm}}{\mathrm{kg}}]",
        color=:red,
        xticks=:none,
        label="x(t)",
        legend=:right
    )
    savefig( "rogers_fig2.svg")


    #wykres 3 i 4
    tspan2 = (0.0, 30.0)
    plt = plot(title="Supersaturation for different initial conditions", xlabel="t", ylabel="S-1 [%]", framestyle = :box)
    plt2 = plot(title="Radius for different initial conditions", xlabel="t", ylabel="Radius "*L" [\mathrm{µm}]", framestyle = :box)


    for s in [1.02,1.01,1.0,0.99, 0.98]
        u0[Int(S)] = s
        prob2 = ODEProblem(rhs!, u0, tspan2,param)
        sol2 = solve(prob2, saveat = 0.1)
        plot!(plt, sol2.t, (sol2[Int(S),:] .- 1)*100, 
        label="S(t=0) = $(round((s-1)*100,digits=2)) %"
        )
        plot!(plt2, sol2.t, sol2[Int(r),:] .* 1e6, 
        label="S(t=0) = $(round((s-1)*100, digits = 2)) %"
        )

    end
    savefig(plt,"rogers_fig3.svg")
    savefig(plt2,"rogers_fig4.svg")


    # FIGURA 5

    plt = Plots.plot(title="Supersaturation for different updraft speed", 
        xlabel="t", 
        ylabel="S-1 [%]"
    )
    u0[Int(S)]=1.0

    for u in [1.0,5.0,15.0,10.0]
        param[Int(U)]= u
        prob2 = ODEProblem(rhs!, u0, tspan,param)
        sol2 = solve(prob2, saveat = 0.1)
        Plots.plot!(plt, sol2.t, (sol2[Int(S),:] .- 1)*100, 
        label="U = $(u) "*L" [\frac{\mathrm{m}}{\mathrm{s}}]"
        )


    end
    savefig(plt,"rogers_fig5.svg")


    # FIGURA6
    plt = Plots.plot(title="Supersaturation for various values \n of droplet concentra-
    tion", 
        xlabel="t", 
        ylabel=LaTeXStrings.L"\mathrm{ln}0(\frac{s_\infty-s}{s_0-s})",
        ylims = (-4,0.5)
    )


    tspan3 = (0.0,6.0) 
    taus = Dict(50.0 => 5.45, 100.0 => 3.17,500.0 => 0.73,200.0 => 1.73)
    colors = Dict(50.0 => "red", 100.0 => "green",500.0 => "blue",200.0 => "deeppink")
    for n in [50.0,100.0,500.0,200.0]
        param[3] =  n*10^(6) / param[2]  # 200 /cm^3 -> 2e6 /m^3
        prob3 = ODEProblem(rhs!, u0, tspan3,param)
        sol3 = solve(prob3, saveat = 0.1)
        s_infty = s_inf.(sol3[Int(T),:],sol3[Int(p),:],sol3[Int(r),:], Ref(param))
        s=sol3[Int(S),:] .- 1
        line = (s_infty.-s)./(s_infty.-s[1])
        safe_line = max.(line, 0)
        line = log.(safe_line)
        Plots.plot!(plt, sol3.t, line, 
        label=L"n_{0} = $(n)\ [\mathrm{cm}^{-3}]",
        color = colors[n]
        )
        exponential = -sol3.t./taus[n]
        Plots.plot!(plt, sol3.t, exponential;
        label = "",         
        linestyle = :dash,   
        color = colors[n]        
        )


    end
    savefig(plt,"rogers_fig6.svg")



end  