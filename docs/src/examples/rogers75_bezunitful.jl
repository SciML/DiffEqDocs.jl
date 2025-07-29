using Plots
using DifferentialEquations
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

# Wektor początkowy u0
u0 = Vector{Float64}(undef, 4)
u0[Int(p)] = 800*100         # 800 mbar = 80000 Pa
u0[Int(T)] = 273.15 + 7.0    # 7°C = 280.15 K
u0[Int(S)] = 1.0             # bezwymiarowe
u0[Int(r)] = 8.0*10^(-6)          # 8 µm = 8e-6 m
param=zeros(2)
param[1]=10.0 # U
param[2]=u0[Int(p)]/C.R/u0[Int(T)]#rho_0
# Zakres czasu 
tspan = (0.0, 20.0)

# Definicja i rozwiązanie problemu ODE
prob = ODEProblem(rhs!, u0, tspan,param)
sol = solve(prob)

# Wykres
plot(sol.t, sol[Int(r),:] .* 1e6, 
    ylabel="Promień [um]", xlabel="Czas [s]"
)
savefig( "test_r.svg")
plot(sol.t, (sol[Int(S),:] .- 1).*100,
    ylabel="Przesycenie [%]",
    xlabel="Czas [s]"
)
savefig( "test_s.svg")

