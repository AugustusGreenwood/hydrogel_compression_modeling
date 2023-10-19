function getNodePositions(feb_file::String)::Dict{Int64, Vector{Float64}}
    id_pos_dict = Dict{Int64, Vector{Float64}}()

    io = open(feb_file, "r")
    while eof(io) == false
        line = strip(readline(io))

        if startswith(line, "<Nodes name=")
            new_object_dict = __parse_feb_file_node_position_section(io)
            merge!(id_pos_dict, new_object_dict)
        end
    end
    return id_pos_dict
end

function getElementConnectivity(feb_file::String)
    id_connectivity_dict = Dict{Int64, Vector{Int64}}()

    io = open(feb_file, "r")
    while eof(io) == false
        line = strip(readline(io))

        if startswith(line, "<Elements type=")
            new_object_dict = __parse_feb_file_element_connectivity_section(io)
            merge!(id_connectivity_dict, new_object_dict)
        end
    end
    return id_connectivity_dict
end

function getElementCenterOfMass(feb_file::String)::Dict{Int64, Vector{Float64}}
    node_positions = getNodePositions(feb_file)
    element_connectivity = getElementConnectivity(feb_file)

    element_COM = Dict{Int64, Vector{Float64}}()
    for (ele_id, conn_vec) in element_connectivity
        order = length(conn_vec)
        
        ele_node_pos_matrix = Matrix{Float64}(undef, 3, order)
        for (i, conn) in enumerate(conn_vec)
            ele_node_pos_matrix[:,i] = node_positions[conn]
        end
        
        element_COM[ele_id] = round.(mean(ele_node_pos_matrix, dims=2)[:], digits=10)
    end
    return element_COM
end

function getData(data_file::String; out_type::Symbol=:dict)
    df = CSV.read(data_file, delim='\t', DataFrame)
    rename!(df, "x" => "x0")

    if out_type == :dataframe 
        return df
    elseif out_type == :dict 
        return __dataframe_to_dict(df)
    else
        display("out_type $out_type not recognized. Returning Dictionary")
        return __dataframe_to_dict(df)
    end
end

function getDataAtIndex(data_file::String, index::Int64)
    data = getData(data_file)

    out = Vector{Float64}(undef, length(data)-1)
    i = 1
    for (node, var) in data
        if node != 0 
            out[i] = var[index]
            i += 1
        end
    end
    return out
end

function getAllDataAtIndex(root::String, index::Int64)
    data_files = data_files = [file for file in readdir("$root/", join=true) if endswith(file, ".dat")]
    data_files = data_files[[7,2,3,4,5,6,1]]
    vec_data = getDataAtIndex.(data_files, index)
    data = reduce(hcat, vec_data)'
    return data
end


function getTimepointData(dir::String)
    timepoints = getTimepointIndices(dir)
    return [getAllDataAtIndex(dir, point) for point in timepoints]
end


function getTimepointIndices(dir::String)
    f = open(joinpath(["$dir", "indices.txt"]))
    readline(f)
    
    r, eor, eor_dwell, eoc, eoc_dwell = @scanf(readline(f), "%i, %i, %i, %i\n", Int64, Int64, Int64, Int64)
    close(f)
    return eor, eor_dwell, eoc, eoc_dwell
end


function __parse_feb_file_node_position_line(line::String)::Tuple{Int64, Vector{Float64}}
    r, id, s = @scanf(line, "\t\t\t<node id=\"%i\">%[^<]s</node>\n", Int64, String)
    pos = parse.(Float64, split(s, ','))
    
    if r != 2
        display("Something went wrong! Expected to read 2 items and $r were actually read.")
    elseif length(pos) != 3
        display("Something went wrong! $pos has incorrect dimentions")
    end

    return id, pos
end

function __parse_feb_file_node_position_section(file::IOStream)
    id_pos_dict = Dict{Int64, Vector{Float64}}()
    while true
        line = readline(file)

        if startswith(line, "\t\t</Nodes>")
            break
        end

        id, pos_vec = __parse_feb_file_node_position_line(line)
        id_pos_dict[id] = pos_vec
    end
    return id_pos_dict
end

function __parse_feb_file_element_connectivity_line(line::String)::Tuple{Int64, Vector{Int64}}
    r, id, string_conn = @scanf(line, "\t\t\t<elem id=\"%i\">%[^<]s</elem>\n", Int64, String)

    if r != 2
        display("Something went wrong! Expected to read 2 items and $r were actually read.")
    end

    return id, parse.(Int64, split(string_conn, ','))
end

function __parse_feb_file_element_connectivity_section(file::IOStream)
    id_connectivity_dict = Dict{Int64, Vector{Int64}}()
    while true
        line = readline(file)

        if startswith(line, "\t\t</Elements>")
            break
        end

        id, id_connectivity_dict[id] = __parse_feb_file_element_connectivity_line(line)
    end
    return id_connectivity_dict
end

function __dataframe_to_dict(dataframe::DataFrame)
    return Dict([parse(Int64, col[2:end]) => dataframe[!, Symbol(col)] for col in names(dataframe)]) 
end

