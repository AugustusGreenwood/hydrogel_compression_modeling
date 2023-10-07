global const EOC_INDEX::Int64, EOC_TIME::Float64, EOR_INDEX::Int64, EOR_TIME::Float64 = __get_eoc_eor_locations()::Tuple{Int64, Float64, Int64, Float64}

function __get_eoc_eor_locations()::Tuple{Int64, Float64, Int64, Float64}
    file = open("locations.dat", "r")
    eoc_line, eor_line = readline(file), readline(file)

    r, eoc_indx, eoc_time = @scanf(eoc_line, "%i,%f", Int64, Float64)
    r, eor_indx, eor_time = @scanf(eor_line, "%i,%f", Int64, Float64)

    return eoc_indx, eoc_time, eor_indx, eor_time
end

function __get_r_h(feb_file::String)
    node_loc = getNodePositions(feb_file)

    node_rh = Dict{Int64, Matrix{Float64}}()
    for (node, loc) in node_loc
        r = sqrt(loc[1]^2 + loc[2]^2)
        h = loc[3]
        node_rh[node] = [r h]
    end
    return node_rh
end


function __get_r_h_data(data_file::String, index::Int64)
    data = getData(data_file; out_type=:dict)
    
    i = 1
    out = Vector{Float64}(undef, length(data)-1)
    for (node, point) in data
        if node != 0
            out[i] = point[index]
            i += 1
        end
    end
    return out
end

function __get_r_h_eoc_eor_data(data_file::String, feb_file::String)
    data = getData(data_file; out_type=:dict)
    node_rh = __get_r_h(feb_file)

    node_rhd = Dict{Int64, Matrix{Float64}}()
    for (node, point) in data
        if node != 0
            node_rhd[node] = [node_rh[node]... point[EOC_INDEX] point[EOR_INDEX]] 
        end
    end
    return node_rhd
end


function __dict_to_mat(data::Dict)
    return reduce(vcat, values(data))
end

function __dict_to_mat_sort(data::Dict)
    return sortslices(reduce(vcat, values(data)), dims=1)
end


function __scaleData(data::AbstractVector{Float64})
    min, max = minimum(data), maximum(data)

    return @. 2 * (data - min) / (max - min) - 1
end

function __scaleData(data::AbstractMatrix{Float64})
    out = similar(data)
    for (i, row) in enumerate(eachrow(data))
        out[i,:] = __scaleData(row)
    end
    return out
end




