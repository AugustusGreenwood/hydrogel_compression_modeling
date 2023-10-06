function __strain_history_summation(time::Float64, betas::Vector{Float64}, sys::System)::Float64
    summation = 0
    for beta in betas
        summation += 4 * exp(-beta^2 * time / sys.tg) / (9*(1 - sys.nu)^2*beta^2 - 8*(1 + sys.nu)*(1 - 2*sys.nu))
    end
    return summation
end

function strain_history(time::Float64, betas::Vector{Float64}, sys::System)::Float64
    -sys.F / (sys.E*sys.A) * (1 - (1 - sys.nu^2)*(1 - 2*sys.nu) * __strain_history_summation(time, betas, sys))
end


function plot_strain_history_validation(data_file::String, sys::System)
    betas = get_characteristic_equation_roots_beta(sys)
    data = get_febio_data_as_matrix(data_file)
    fe_time, fe_strain = data[:,1], -data[:,2]/sys.h
    analytical = strain_history.(fe_time, Ref(betas), sys)


    plot(fe_time, analytical)
    scatter!(fe_time, fe_strain)
end

function plot_strain_history_validation!(data_file::String, sys::System)
    betas = get_characteristic_equation_roots_beta(sys)
    data = get_febio_data_as_matrix(data_file)
    fe_time, fe_strain = data[:,1], -data[:,2]/sys.h
    analytical = strain_history.(fe_time, Ref(betas), sys)

    plot!(fe_time, analytical)
    scatter!(fe_time, fe_strain)
end



function output_strain_history_validation(name::String, data_file::String, sys::System)
    betas = get_characteristic_equation_roots_beta(sys)
    data = get_febio_data_as_matrix(data_file)
    fe_time, fe_strain = data[:,1], -data[:,2]/sys.h
    analytical = strain_history.(fe_time, Ref(betas), sys)


    writedlm("$name", [fe_time/sys.tg analytical fe_strain], ',')
end