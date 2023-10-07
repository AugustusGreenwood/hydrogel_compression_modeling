function outputVectorData(data_file::String, feb_file::String; seperate::Bool=false)
    root, out = splitdir(data_file)
    rhd_dict = __get_r_h_eoc_eor_data(data_file, feb_file)
	rhd_mat = __dict_to_mat_sort(rhd_dict)
    if seperate == false
        writedlm("$root/r_h_eoc_eor_$out", rhd_mat', ',')
    else
        r_eoc_mat = [rhd_mat[:,1] rhd_mat[:,3]]
        h_eoc_mat = [rhd_mat[:,2] rhd_mat[:,3]]
        r_eor_mat = [rhd_mat[:,1] rhd_mat[:,4]]
        h_eor_mat = [rhd_mat[:,2] rhd_mat[:,4]]

        writedlm("$root/r_eoc_$out", r_eoc_mat, ',')
        writedlm("$root/h_eoc_$out", h_eoc_mat, ',')
        writedlm("$root/r_eor_$out", r_eor_mat, ',')
        writedlm("$root/h_eor_$out", h_eor_mat, ',')
    end
end



function outputMatrixData(data_file::String, feb_file::String)
    root, out = splitdir(data_file)
    rhd_dict = __get_r_h_eoc_eor_data(data_file, feb_file)
    rhd_mat = __dict_to_mat_sort(rhd_dict)
    r = unique(rhd_mat[:,1])
    h = unique(rhd_mat[:,2])
    eoc = reshape(rhd_mat[:,3], 26, 26)
    eor = reshape(rhd_mat[:,4], 26, 26)
    eoc_mat = hcat([0;h], vcat(r', eoc))
    eor_mat = hcat([0;h], vcat(r', eor))
    writedlm("$root/eoc_$out", eoc_mat, ',')
    writedlm("$root/eor_$out", eor_mat, ',')
end


function outputAllData(dir::String, feb_file::String)
    in_files = [joinpath(dir, file) for file in readdir(dir) if endswith(file, ".dat")]
    Threads.@threads for file in in_files
        outputVectorData(file, feb_file; seperate=true)
        outputMatrixData(file, feb_file)
    end
end


function scaleData(data::AbstractVector{Float64})
    min, max = minimum(data), maximum(data)

    return @. 2 * (data - min) / (max - min) - 1
end

function scaleData(data::AbstractMatrix{Float64})
    out = similar(data)
    for (i, row) in enumerate(eachrow(data))
        out[i,:] = scaleData(row)
    end
    return out
end





function __dict_to_mat(data::Dict)
    return reduce(vcat, values(data))
end

function __dict_to_mat_sort(data::Dict)
    return sortslices(reduce(vcat, values(data)), dims=1)
end


