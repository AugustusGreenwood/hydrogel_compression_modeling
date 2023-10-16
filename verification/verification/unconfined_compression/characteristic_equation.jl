function characteristic_equation_alpha(x::Number, nu::Float64)::Float64
    return besselj1(x) - (1 - nu) * x * besselj0(x) / (1 - 2*nu)
end

function __get_characteristic_equation_roots_alpha(start::Number, stop::Number, nu::Float64)::Vector{Float64}
    f(x) = characteristic_equation_alpha(x, nu)
    return find_zeros(f, start, stop, no_pts=Int(round((2*stop-2*start), sigdigits=1)))
end


function get_characteristic_equation_roots_alpha(sys::System)::Vector{Float64}
    alphas = CSV.read("characteristic_equation_roots_alpha.CSV", DataFrame)
    return alphas[!, "$(sys.nu)"]
end

function get_characteristic_equation_roots_alpha(nu::Float64)::Vector{Float64}
    alphas = CSV.read("characteristic_equation_roots_alpha.CSV", DataFrame)
    return alphas[!, "$(nu)"]
end




function characteristic_equation_beta(x::Number, nu::Float64)::Float64
    return besselj0(x) - 4 * (1 - 2*nu) * besselj1(x) / (3 * (1 - nu) * x)
end

function __get_characteristic_equation_roots_beta(start::Number, stop::Number, nu::Float64)::Vector{Float64}
    f(x) = characteristic_equation_beta(x, nu)
    return find_zeros(f, start, stop, no_pts=Int(round((2*stop-2*start), sigdigits=1)))
end


function get_characteristic_equation_roots_beta(sys::System)::Vector{Float64}
    alphas = CSV.read("characteristic_equation_roots_beta.CSV", DataFrame)
    return alphas[!, "$(sys.nu)"]
end

function get_characteristic_equation_roots_beta(nu::Float64)::Vector{Float64}
    alphas = CSV.read("characteristic_equation_roots_beta.CSV", DataFrame)
    return alphas[!, "$(nu)"]
end