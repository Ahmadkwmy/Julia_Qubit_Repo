using Plots
using DifferentialEquations
using OrdinaryDiffEq

# Define the system equations
function TwoLevelSystem(dρ::Vector{ComplexF64}, ρ::Vector{ComplexF64}, params, t)
    Ω, Δ, γ, Γ = params
    dρ[1] =  Γ * ρ[4] + (im * Ω) / 2 * (ρ[2] - ρ[3])
    dρ[4] = -Γ * ρ[4] - (im * Ω) / 2 * (ρ[2] - ρ[3])
    dρ[2] = -(γ + im * Δ) * ρ[2] - (im * Ω) / 2 * (ρ[4] - ρ[1])
    dρ[3] = -(γ - im * Δ) * ρ[3] + (im * Ω) / 2 * (ρ[4] - ρ[1])
end


ρ_0 = ComplexF64[1.0 + 0im, 0.0 + 0im, 0.0 + 0im, 0.0 + 0im]

params = (1.0, 0.0, 0.1, 0.1)

time_range = (0.0, 10.0)

probODE = ODEProblem(TwoLevelSystem, ρ_0, time_range, params)

sol = solve(probODE, Tsit5())

t = sol.t
ρ11 = real.(sol[1, :])
ρ22 = real.(sol[4, :])

plot(t, ρ11)
plot!(t, ρ22)
