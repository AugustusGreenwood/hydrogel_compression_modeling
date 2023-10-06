function get_pressure_summation(time::Float64, depth::Float64, sys::System)::Float64
    summation = 0
    for n in 1:1000
        λ = calculate_lambda(n)
        summation += (1 / λ) * sin(λ * depth / sys.h) * exp(-λ^2 * time / sys.tg)
    end
    return summation
end


function calculate_pressure(time::Float64, depth::Float64, sys::System)
    return 2 * sys.F / sys.A * get_pressure_summation(time, depth, sys)
end


function plot_non_dim_pressure_vary_depth(sys::System, time::Float64; plot_options...)
    depth = LinRange(0, sys.h, 1000)
    p = calculate_pressure.(time, depth, sys)
    normalized_p = p * sys.A / sys.F

    plot(depth/sys.h, normalized_p; plot_options...)
end

function plot_non_dim_pressure_vary_time(sys::System, depth::Float64; plot_options...)
    time = LinRange(0, sys.tg*5, 100000)
    p = calculate_pressure.(time, depth, sys)
    normalized_p = p * sys.A / sys.F

    plot(time/sys.tg, normalized_p; plot_options...)
end

function plot_non_dim_pressure_vary_depth!(sys::System, time::Float64; plot_options...)
    depth = LinRange(0, sys.h, 1000)
    p = calculate_pressure.(time, depth, sys)
    normalized_p = p * sys.A / sys.F

    plot!(depth/sys.h, normalized_p; plot_options...)
end

function plot_non_dim_pressure_vary_time!(sys::System, depth::Float64; plot_options...)
    time = LinRange(0, sys.tg*5, 100000)
    p = calculate_pressure.(time, depth, sys)
    normalized_p = p * sys.A / sys.F

    plot!(time/sys.tg, normalized_p; plot_options...)
end


function plot_time_pressure_validation(data_file::String, feb_file::String, sys::System)
    data = get_febio_data(data_file)
    node_pos = get_node_positions(feb_file)
    node_depths = get_node_depths(sys, node_pos, data)
    time_norm = data["x"] / sys.tg

    plot(size=(800, 700))
    for (i, (k, v)) in enumerate(node_depths)
        plot_non_dim_pressure_vary_time!(sys, v; label="Depth $v")
        p_norm = data[k] * sys.A / sys.F

        scatter!(time_norm[2:5:end], -p_norm[2:5:end], label="Depth $v", linestyle=:dot, markersize=6, markerstrokealpha=1, markeralpha=0.5)
    end
    plot!()
end


function plot_height_pressure_validation(index::Int64, data_file::String, feb_file::String, sys::System; plot_options...)
    data = get_febio_data(data_file)
    node_pos = get_node_positions(feb_file)
    node_depths = get_node_depths(sys, node_pos, data)
    time = data["x"][index]

    plot(size=(800, 700), legend=false)
    heights = Vector{Float64}(undef, length(node_depths))
    pressures = Vector{Float64}(undef, length(node_depths))
    for (i, (k, v)) in enumerate(node_depths)
        heights[i] = v / sys.h
        pressures[i] = -data[k][index] * sys.A / sys.F
    end

    plot_non_dim_pressure_vary_depth(sys, time; label="Analytical time $time", plot_options...)    

    scatter!(heights, pressures, label="FEM time $time", linestyle=:dot, markersize=6, markerstrokealpha=1, markeralpha=0.5; plot_options...)
end


function plot_height_pressure_validation!(index::Int64, data_file::String, feb_file::String, sys::System; plot_options...)
    data = get_febio_data(data_file)
    node_pos = get_node_positions(feb_file)
    node_depths = get_node_depths(sys, node_pos, data)
    time = data["x"][index]

    heights = Vector{Float64}(undef, length(node_depths))
    pressures = Vector{Float64}(undef, length(node_depths))
    for (i, (k, v)) in enumerate(node_depths)
        heights[i] = v / sys.h
        pressures[i] = -data[k][index] * sys.A / sys.F
    end

    plot_non_dim_pressure_vary_depth!(sys, time; label="Analytical time $time", plot_options...)    

    scatter!(heights, pressures, label="FEM time $time", linestyle=:dot, markersize=6, markerstrokealpha=1, markeralpha=0.5; plot_options...)
end










function output_non_dim_pressure_vary_depth(sys::System, time::Float64)
    depth = LinRange(0, sys.h, 1000)
    p = calculate_pressure.(time, depth, sys)
    normalized_p = p * sys.A / sys.F

    writedlm("theory_non_dim_pressure_vary_depth__time$time.csv", [depth normalized_p], ',')
end

function output_non_dim_pressure_vary_time(sys::System, depth::Float64)
    time = LinRange(0, sys.tg*5, 100000)
    p = calculate_pressure.(time, depth, sys)
    normalized_p = p * sys.A / sys.F

    writedlm("thoery_non_dim_pressure_vary_time__depth$depth.csv", [time normalized_p], ',')
end




function output_height_pressure_validation(index::Int64, data_file::String, feb_file::String, sys::System; plot_options...)
    data = get_febio_data(data_file)
    node_pos = get_node_positions(feb_file)
    node_depths = get_node_depths(sys, node_pos, data)
    time = data["x"][index]


    heights = Vector{Float64}(undef, length(node_depths))
    pressures = Vector{Float64}(undef, length(node_depths))
    for (i, (k, v)) in enumerate(node_depths)
        heights[i] = v / sys.h
        pressures[i] = data[k][index] * sys.A / sys.F
    end

    p = calculate_pressure.(time, heights, sys)
    normalized_p = p * sys.A / sys.F

    writedlm("non_dim_pressure_vary_depth__time$time.csv", [heights normalized_p -pressures], ',')
end

function output_height_pressure_validation(indices::StepRange{Int64, Int64}, data_file::String, feb_file::String, sys::System; plot_options...)
    data = get_febio_data(data_file)
    node_pos = get_node_positions(feb_file)
    node_depths = get_node_depths(sys, node_pos, data)
    outdata = ["Normalized height"]
    outdata = [outdata; [v / sys.h for (k,v) in node_depths]]

    for index in indices
        time = data["x"][index]

        heights = Vector{Float64}(undef, length(node_depths))
        pressures = Vector{Float64}(undef, length(node_depths))
        for (i, (k, v)) in enumerate(node_depths)
            heights[i] = v / sys.h
            pressures[i] = data[k][index] * sys.A / sys.F
        end

        p = calculate_pressure.(time, heights, sys)
        normalized_p = p * sys.A / sys.F

        theory = ["Theory time $(time)s"; normalized_p]
        fem = ["FEM time $(time)s"; -pressures]

        outdata = [outdata theory]
        outdata = [outdata fem]
    end

    writedlm("non_dim_pressure_vary_depth", outdata, ',')
end





function output_time_pressure_validation(data_file::String, feb_file::String, sys::System)
    data = get_febio_data(data_file)
    node_pos = get_node_positions(feb_file)
    node_depths = get_node_depths(sys, node_pos, data)
    time_norm = data["x"] / sys.tg
    outdata = ["Nondimentional time"; time_norm]

    for (k, v) in node_depths
        p = calculate_pressure.(time_norm * sys.tg, v, sys)
        normalized_p = p * sys.A / sys.F
        p_norm = data[k] * sys.A / sys.F; 

        theory_out = ["Theory depth $(v)mm"; normalized_p]
        fem_out = ["FEM depth $(v)mm"; -p_norm]

        outdata = [outdata theory_out fem_out]    
    end

    writedlm("non_dim_pressure_vary_time.csv", outdata, ',')
end
