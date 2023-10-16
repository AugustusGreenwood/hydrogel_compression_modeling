using Roots
using Plots; plotlyjs()
using SpecialFunctions
using DataFrames, CSV
using DelimitedFiles
using Distributed

struct System
    E::Float64
    nu::Float64
    k::Float64
    ϕs::Float64
    r::Float64
    h::Float64
    θ::Float64
    ϵ::Float64
    F::Float64
    A::Float64
    Ha::Float64
    G::Float64
    tg::Float64
    ϕf::Float64
    alpha::Float64
    function System(E::Float64, nu::Float64, k::Float64, ϕs::Float64, r::Float64, h::Float64, θ::Float64, ϵ::Float64, F::Float64)
        A = π * r^2 / (360 / θ)
        Ha = calculate_Ha(E, nu)
        G = calculate_G(E, nu)
        tg = calculate_tg(E, nu, r, k)
        ϕf = 1 - ϕs
        alpha = ϕs / ϕf
        return new(E, nu, k, ϕs, r, h, θ, ϵ, F, A, Ha, G, tg, ϕf, alpha)
    end
    function System(E::Float64, nu::Float64, k::Float64, ϕs::Float64, r::Float64, h::Float64, θ::Float64, ϵ::Float64, F::Float64, slices::Int64)
        A = sin(deg2rad(θ/slices)) * r^2 * slices / 2
        Ha = calculate_Ha(E, nu)
        G = calculate_G(E, nu)
        tg = calculate_tg(E, nu, r, k)
        ϕf = 1 - ϕs
        alpha = ϕs / ϕf
        return new(E, nu, k, ϕs, r, h, θ, ϵ, F, A, Ha, G, tg, ϕf, alpha)
    end
end

Base.Broadcast.broadcastable(sys::System) = Ref(sys)

function calculate_Ha(E::Float64, nu::Float64)::Float64
    return E * (1 - nu) / ((1 + nu) * (1 - 2*nu))
end

function calculate_G(E::Float64, nu::Float64)::Float64
    return E / (2 * (1 + nu))
end

function calculate_tg(E::Float64, nu::Float64, r::Float64, k::Float64)::Float64
    Ha = calculate_Ha(E, nu)
    return r^2 / (Ha * k)
end

include("characteristic_equation.jl")
include("radial_displacement.jl")
include("fluid_pressure.jl")
include("strain_history.jl")
include("../FEBioInteract.jl")






