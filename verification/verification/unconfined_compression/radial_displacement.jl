function __prescribed_strain_radial_displacement_summation(time::Float64, alphas::Vector{Float64}, sys::System)::Float64
    summation = 0
    for alpha in alphas
        summation += exp(-alpha^2 * time / sys.tg) / (alpha^2*(1 - sys.nu)^2 - (1 - 2*sys.nu))
    end
    return summation
end

function __prescribed_strain_radial_displacement(time::Float64, alphas::Vector{Float64}, sys::System)::Float64
    return sys.r * sys.ϵ * (sys.nu + (1 - 2*sys.nu)*(1 - sys.nu) * __prescribed_strain_radial_displacement_summation(time, alphas, sys))
end


function __prescribed_force_radial_displacement_summation(time::Float64, betas::Vector{Float64}, sys::System)::Float64
    summation = 0
    for beta in betas
        summation += 4 * exp(-beta^2 * time / sys.tg) / (9*(1 - sys.nu)^2*beta^2 - 8*(1 + sys.nu)*(1 - 2*sys.nu))
    end
    return summation
end

function __prescribed_force_radial_displacement(time::Float64, betas::Vector{Float64}, sys::System)::Float64
    return -sys.F * sys.r / sys.A * (sys.nu/sys.E + (1 - sys.nu^2)*(1 - 2*sys.nu) * __prescribed_force_radial_displacement_summation(time, betas, sys) / sys.E)
end

function radial_displacement(sys::System)::Vector{Float64}
    time = LinRange(0, 2*sys.tg, 100000)

    if sys.ϵ == 0
        betas = get_characteristic_equation_roots_beta(sys)
        return __prescribed_force_radial_displacement.(time, Ref(betas), sys)
    else
        alphas = get_characteristic_equation_roots_alpha(sys)
        return __prescribed_strain_radial_displacement.(time, Ref(alphas), sys)
    end
end


function radial_displacement(time::Float64, sys::System)::Float64
    if sys.ϵ == 0
        betas = get_characteristic_equation_roots_beta(sys)
        return __prescribed_force_radial_displacement(time, betas, sys)
    else
        alphas = get_characteristic_equation_roots_alpha(sys)
        return __prescribed_strain_radial_displacement(time, alphas, sys)
    end
end

function plot_radial_displacement_validation(data_file::String, sys::System)
    data = get_febio_data_as_matrix(data_file)
    fe_time, fe_strain = data[:,1], data[:,2] 
    analytical = radial_displacement.(fe_time, sys)


    plot(fe_time/sys.tg, analytical)
    scatter!(fe_time/sys.tg, fe_strain)
end


function plot_radial_displacement_validation!(data_file::String, sys::System)
    data = get_febio_data_as_matrix(data_file)
    fe_time, fe_strain = data[:,1], data[:,2] 
    analytical = radial_displacement.(fe_time, sys)

    plot!(fe_time/sys.tg, analytical)
    scatter!(fe_time/sys.tg, fe_strain)
end


function output_radial_displacement_validation(name::String, data_file::String, sys::System)
    data = get_febio_data_as_matrix(data_file)
    fe_time, fe_strain = data[:,1], data[:,2] 
    analytical = radial_displacement.(fe_time, sys)

    writedlm("$name", [fe_time/sys.tg fe_strain analytical], ',')
end

