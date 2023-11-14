function readFEBioDataFile(path::String)::Dict{Int64,Vector{Float64}}
    matrix_data, keys = open(path, "r") do file
        return _get_FEBio_data_and_keys(file)
    end
    # Wanted to speed this up, and concurrent write to a dict is a no-no. So, 
    # each loop we lock data_dict, write and unclock after
    lk = ReentrantLock()
    data_dict = Dict{Int64,Vector{Float64}}()
    Threads.@threads for (i, key) in collect(enumerate(keys))
        lock(lk) do
            data_dict[key] = matrix_data[i, :]
        end
    end
    return data_dict
end

#=
Loads all the data at what ever time or index is given into a vector. Might be
a better way to do getting the data_files, but not sure
=#
function getAllFEBioDataAtTime(root_dir::String, time::Union{Float64,Int64})
    data_file_names = [
        "hpos.dat", "rpos.dat", "lagrange_strain_invariant_1.dat",
        "lagrange_strain_invariant_2.dat", "lagrange_strain_invariant_3.dat",
        "fluid_flux_magnitude.dat", "hydrostatic_pressure.dat"
    ]
    paths = joinpath.(root_dir, data_file_names)
    return getFEBioDataAtTime.(paths, time)
end

#=
Maybe this isn't the best way to go about this, but this converts the 
Vector{Vector{Number}} data into a matrix to feed into algorithsm. It may be 
advatangeous to have the original nested vector so I didn't what to just put 
this in getAllFEBioDataAtTime so its seperate. 

I also will almost always use it with a certain data as rows, but may as well 
give an option to change it.
=#
function collectFEBioData(data::Vector{Vector{Float64}}, as_row::Bool=true)::Matrix{Float64}
    matrix_data = reduce(hcat, data)
    if as_row == true
        return matrix_data'
    else
        return matrix_data
    end
end

#=
I don't actually know if I'll ever use this, but it was easy to implement so
=#
function getFEBioDataAtLocation(path::String, index::Int64)::Matrix{Float64}
    matrix_data, _ = open(path, "r") do file
        return _get_FEBio_data_and_keys(file)
    end

    return matrix_data[index, :]
end

#=
These can get a slice of the febio data at some time. This can be done at the
number index, or a specific time. There might be a better way to do this than 
multiple dispatch based whether argumnet two is an int or float but for now 
this is what it is
=#
function getFEBioDataAtTime(path::String, index::Int64)::Vector{Float64}
    matrix_data, _ = open(path, "r") do file
        return _get_FEBio_data_and_keys(file)
    end

    return matrix_data[2:end, index]
end

function getFEBioDataAtTime(path::String, time::Float64)::Vector{Float64}
    matrix_data, _ = open(path, "r") do file
        return _get_FEBio_data_and_keys(file)
    end

    index = findfirst(x -> x == time, matrix_data[1, :])
    if isnothing(index)
        throw(ArgumentError("Given time ($time) is not in data"))
    end
    return matrix_data[2:end, index]
end

#=
Spilts a whole FEBio data file into the keys (first line) and the data 
(everything else) and returns a matrix of that data and the keys. Used to make
it easy for convert the data into a dictionary
=#
function _get_FEBio_data_and_keys(file::IOStream)::Tuple{Matrix{Float64},Vector{Int64}}
    # Read whole file, remove any empty lines, break apart lines by tab character
    nested_vector_string_data = split.(filter(!isempty, readlines(file)), '\t')
    # Get keys (first line), convert x to 0, parse into integers for dictionary
    keys = [key[2:end] for key in popfirst!(nested_vector_string_data)]
    keys[keys.==""] .= "0"
    keys = parse.(Int64, keys)
    # Parse each line from vector of strings to vector of floats
    nested_float_vector_data = [parse.(Float64, vec) for vec in nested_vector_string_data]
    # Convert the Vector{Vector{Float64}} into a matrix
    matrix_float_data = reduce(hcat, nested_float_vector_data)
    return matrix_float_data, keys
end



