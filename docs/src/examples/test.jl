using DifferentialEquations
using Unitful
using Base
using Test
using NamedTupleTools
using Plots

functions = (
    rho = (constants, pressure, temperature) -> pressure / constants.specific_gas_constant / temperature,
    calculate_fk_fd = (constants, temperature, pressure) -> begin
        @show Unitful.dimension(temperature)
        @show Unitful.dimension(constants.convection_constant)
        @show Unitful.dimension((393.0*u"K" / (temperature + 120.0*u"K")) * (temperature / 273.0*u"K")^(3/2))
        
        D = constants.diffusion_consant * temperature / pressure * (393.0*u"K" / (temperature + 120.0*u"K")) * (temperature / (273.0*u"K"))^(3/2)
        # @show Unitful.dimension(D)
        K = constants.convection_constant * (393.0*u"K" / (temperature + 120.0*u"K")) * (temperature / (273.0*u"K"))^(3/2)
        @show Unitful.dimension(K)
        Fk = constants.L^2 * constants.epsilon * constants.rho_l / (K * constants.specific_gas_constant * temperature^2)
        @show Unitful.dimension(Fk)
        Fd = constants.specific_gas_constant * temperature * constants.rho_l / (constants.epsilon * D * es(temperature))
        # @show Unitful.dimension(Fd)
        (Fk, Fd)
    end,
    calculate_q1_q2 = (constants,temperature, pressure) -> begin
      
      Q1 = constants.L * constants.g * constants.epsilon / (constants.specific_gas_constant * constants.cp * temperature^2) - constants.g / (constants.specific_gas_constant * temperature)
      Q2 = constants.specific_gas_constant * temperature/ (constants.epsilon * es(temperature)) + constants.epsilon * constants.L^2 / (constants.cp * temperature * pressure)
      (Q1,Q2)
    end,
    G = (constants,functions, temperature, pressure) -> begin
        Q1,Q2 = functions.calculate_q1_q2(temperature,pressure)
        # @show Unitful.dimension(Q1)
        # @show Unitful.dimension(Q2)
        Fk, Fd  = functions.calculate_fk_fd(temperature,pressure)
        # @show Unitful.dimension(Fk)
        # @show Unitful.dimension(Fd)
        Q1/Q2*(Fk+Fd)/(4*pi*constants.rho_l*functions.rho(pressure,temperature))
    end
    )


constants_u = (
    specific_gas_constant = 287.053 * u"J/(kg*K)",
    diffusion_consant =  uconvert(u"Pa*m^2/(K*s)", 8.28*10^2 * u"dyn/(K*s)"),
    g = 9.81 * u"m/(s^2)",
    epsilon = 0.622 * u"1",
    rho_l = 1000.0 * u"kg/(m^3)",
    cp = 1005 * u"J/(kg*K)",
    L = 2.5e6 * u"J/kg",
    convection_constant = uconvert(u"W/(m*K)", 2.42*10^3 * u"erg/(cm*s*K)")

)
function es(T)
    e = 2.75 * 10^(12) * exp(-(5.44 * 10^3 *u"K") / T)
    return e*u"mbar"
end

function sigma(S, T, p)
    Fk, Fd = calculate_fk_fd(T,p)
    return (S - 1.0) / (Fk + Fd)
end
function s_inf(T,p, r,param)
    G_p_T = G(T,p)
    return G_p_T*param[1]/r/param[3]
end

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

    Z = [uconvert(u"s/g",functions_u.G(functions_u,T,p)) for p in p_values, T in T_values]  # wiersze: T, kolumny: p
    plt =Plots.contour( T_values,p_values, Z, 
            xlabel="T", ylabel="p", yflip=true,
            title="Izolinie funkcji G(p,T)")
    savefig(plt,"rogers_fig7.svg")
end

begin
    constants_n = map(ustrip, constants_u)
    functions_n = mapvalues(f -> Base.Fix1(f, constants_n), functions)

#     solver(functions, constants)
end  
