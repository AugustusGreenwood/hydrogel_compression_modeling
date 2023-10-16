function __prescribed_strain_fluid_pressure(r_hat::Float64, time::Float64, sys::System)::Float64
    if r_hat == 1
        return 0
    else
        error_func_value = 1/sqrt(r_hat) * erfc((1 - r_hat)/2 * sqrt(sys.tg/time))
    end
    return sys.E * sys.ϵ * (1 + sys.alpha) / (2*(1 + sys.nu)) * (1 - error_func_value)
end



function __prescribed_force_fluid_pressure_summation(r_hat::Float64, time::Float64, betas::Vector{Float64}, sys::System)::Float64
    numerator(beta) = 1 - besselj0(beta*r_hat) / besselj0(beta)
    denominator(beta) = (9/8)*(1 - sys.nu)^2*beta^2 - (1 + sys.nu)*(1 - 2*sys.nu)

    summation = 0
    for beta in betas
        summation += numerator(beta) * exp(-beta^2 * time / sys.tg) / denominator(beta)
    end
    return summation
end

function __prescribed_force_fluid_pressure(r_hat::Float64, time::Float64, betas::Vector{Float64}, sys::System)::Float64
    return sys.F * (1 - sys.nu)*(1 - 2*sys.nu) / (sys.A * (1 + sys.alpha)) * __prescribed_force_fluid_pressure_summation(r_hat, time, betas, sys)
end


function plot_time_pressure_validation(data_file::String, feb_file::String, sys::System)
    data = get_febio_data(data_file)
    node_pos = get_node_positions(feb_file)
    node_depths = get_node_radius(node_pos, data)
    time_norm = data["x"] / sys.tg
    betas = get_characteristic_equation_roots_beta(sys)
    #colors = [:red, :black, :blue, :green, :cyan, :purple, :grey, :tan, :pink, :yellow, :brown, :orange]


    plot(size=(1200, 800))
    for (i, (k, v)) in enumerate(node_depths)
        r_hat = v/sys.r
        p = __prescribed_force_fluid_pressure.(r_hat, data["x"], Ref(betas), sys)
        plot!(time_norm, p*(1+sys.alpha), label="Depth $r_hat")

        scatter!(time_norm[2:5:end], data[k][2:5:end], label="Depth $r_hat", linestyle=:dot, markersize=6, markerstrokealpha=1, markeralpha=0.1)
    end
    plot!()
end

function plot_r_hat_pressure_validation(time::Float64, data_file::String, feb_file::String, sys::System)
    data = get_febio_data(data_file)
    node_pos = get_node_positions(feb_file)
    node_depths = get_node_radius(node_pos, data)
    index = findfirst(x->x==time, data["x"])
    betas = get_characteristic_equation_roots_beta(sys)

    plot(size=(1000, 600), legend=false)
    r_hats = Vector{Float64}(undef, length(node_depths))
    pressures = Vector{Float64}(undef, length(node_depths))
    for (i, (k, v)) in enumerate(node_depths)
        r_hats[i] = v/sys.r
        pressures[i] = data[k][index]
    end

    p = __prescribed_force_fluid_pressure.(0:0.001:1, time, Ref(betas), sys)
    plot!(0:0.001:1, p*(1+sys.alpha), label="Time $time")

    scatter!(r_hats, pressures, label="Time $time", linestyle=:dot, markersize=6, markerstrokealpha=1, markeralpha=0.1)
    plot!()
end



function plot_r_hat_pressure_validation(index::Int64, data_file::String, feb_file::String, sys::System)
    data = get_febio_data(data_file)
    node_pos = get_node_positions(feb_file)
    node_depths = get_node_radius(node_pos, data)
    betas = get_characteristic_equation_roots_beta(sys)
    time = data["x"][index]

    plot(size=(800, 600), legend=false)
    r_hats = Vector{Float64}(undef, length(node_depths))
    pressures = Vector{Float64}(undef, length(node_depths))
    for (i, (k, v)) in enumerate(node_depths)
        r_hats[i] = v/sys.r
        pressures[i] = data[k][index]
    end

    p = __prescribed_force_fluid_pressure.(0:0.01:1, time, Ref(betas), sys)
    plot!(0:0.01:1, p*(1+sys.alpha), label="Time $time")

    scatter!(r_hats, pressures, label="Time $time", linestyle=:dot, markersize=6, markerstrokealpha=1, markeralpha=0.1)
end



function plot_time_pressure_validation!(data_file::String, feb_file::String, sys::System)
    data = get_febio_data(data_file)
    node_pos = get_node_positions(feb_file)
    node_depths = get_node_radius(node_pos, data)
    time_norm = data["x"] / sys.tg
    betas = get_characteristic_equation_roots_beta(sys)
    colors = [:red, :black, :blue, :green, :cyan, :purple, :grey, :tan, :pink, :yellow, :brown, :orange]

    for (i, (k, v)) in enumerate(node_depths)
        r_hat = v/sys.r
        p = __prescribed_force_fluid_pressure.(r_hat, data["x"], Ref(betas), sys)
        plot!(time_norm, p*(1+sys.alpha), label="Depth $r_hat", color=colors[i])

        scatter!(time_norm[2:5:end], data[k][2:5:end], label="Depth $r_hat", linestyle=:dot, color=colors[i], markersize=6, markerstrokealpha=1, markeralpha=0.1)
    end
end

function plot_r_hat_pressure_validation!(time::Float64, data_file::String, feb_file::String, sys::System)
    data = get_febio_data(data_file)
    node_pos = get_node_positions(feb_file)
    node_depths = get_node_radius(node_pos, data)
    index = findfirst(x->x==time, data["x"])
    betas = get_characteristic_equation_roots_beta(sys)

    r_hats = Vector{Float64}(undef, length(node_depths))
    pressures = Vector{Float64}(undef, length(node_depths))
    for (i, (k, v)) in enumerate(node_depths)
        r_hats[i] = v/sys.r
        pressures[i] = data[k][index]
    end

    p = __prescribed_force_fluid_pressure.(sort(r_hats), time, Ref(betas), sys)
    plot!(sort(r_hats), p*(1+sys.alpha), label="Time $time")
    scatter!(r_hats, pressures, label="Depth $time", linestyle=:dot, markersize=6, markerstrokealpha=1, markeralpha=0.1)
end


function plot_r_hat_pressure_validation!(index::Int64, data_file::String, feb_file::String, sys::System)
    data = get_febio_data(data_file)
    node_pos = get_node_positions(feb_file)
    node_depths = get_node_radius(node_pos, data)
    betas = get_characteristic_equation_roots_beta(sys)
    time = data["x"][index]

    r_hats = Vector{Float64}(undef, length(node_depths))
    pressures = Vector{Float64}(undef, length(node_depths))
    for (i, (k, v)) in enumerate(node_depths)
        r_hats[i] = v/sys.r
        pressures[i] = data[k][index]
    end

    p = __prescribed_force_fluid_pressure.(0:0.01:1, time, Ref(betas), sys)
    plot!(0:0.01:1, p*(1+sys.alpha), label="Time $time")

    scatter!(r_hats, pressures, label="Time $time", linestyle=:dot, markersize=6, markerstrokealpha=1, markeralpha=0.1)
    plot!()
end





function output_r_hat_pressure_validation(index::Int64, data_file::String, feb_file::String, sys::System)
    data = get_febio_data(data_file)
    node_pos = get_node_positions(feb_file)
    node_depths = get_node_radius(node_pos, data)
    betas = get_characteristic_equation_roots_beta(sys)
    time = data["x"][index]


    r_hats = Vector{Float64}(undef, length(node_depths))
    pressures = Vector{Float64}(undef, length(node_depths))
    for (i, (k, v)) in enumerate(node_depths)
        r_hats[i] = v/sys.r
        pressures[i] = data[k][index]
    end

    p = __prescribed_force_fluid_pressure.(r_hats, time, Ref(betas), sys) * (1 + sys.alpha)
    writedlm("r_hat_fluid_presssure__time$time.CSV", [r_hats p pressures], ',')
end


function output_r_hat_pressure_validation(indices::StepRange{Int64, Int64}, data_file::String, feb_file::String, sys::System)
    data = get_febio_data(data_file)
    node_pos = get_node_positions(feb_file)
    node_depths = get_node_radius(node_pos, data)
    betas = get_characteristic_equation_roots_beta(sys)
    outdata = ["Normalized radial position"]
    outdata = [outdata; [v/sys.r for (k, v) in node_depths]]


    for index in indices
        time = data["x"][index]

        r_hats = Vector{Float64}(undef, length(node_depths))
        pressures = Vector{Float64}(undef, length(node_depths))
        for (i, (k, v)) in enumerate(node_depths)
            r_hats[i] = v/sys.r
            pressures[i] = data[k][index]
        end

        p = __prescribed_force_fluid_pressure.(r_hats, time, Ref(betas), sys) * (1 + sys.alpha)

        theory = ["Theory $(time)s"; p]
        fem = ["FEM $(time)s"; pressures]

        outdata = [outdata theory fem]
    end
    writedlm("r_hat_fluid_presssure__time$time.CSV", outdata, ',')
end




function output_time_pressure_validation(data_file::String, feb_file::String, sys::System)
    data = get_febio_data(data_file)
    node_pos = get_node_positions(feb_file)
    node_depths = get_node_radius(node_pos, data)
    time_norm = data["x"] / sys.tg
    betas = get_characteristic_equation_roots_beta(sys)
    outdata = ["Nondimentionalized time"; time_norm]


    for (k, v) in node_depths
        r_hat = v/sys.r
        p = __prescribed_force_fluid_pressure.(r_hat, data["x"], Ref(betas), sys)*(1+sys.alpha)

        theory = ["Theory radial position $(v)mm"; p]
        fem = ["FEM radial position $(v)mm"; data[k]]

        outdata = [outdata theory fem]
    end

    writedlm("time_fluid_pressure.CSV", outdata, ',')
end


