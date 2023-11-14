using Scanf

function get_febio_data(file::String)::Dict{String, Vector{Float64}}
    io = open(file, "r")
    
    rows = countlines(io) - 1; seek(io, 0)
    
    data = Dict{String, Vector{Float64}}()
    index = Dict{Int64, String}()
    firstline = split(readline(io), '\t')
    for (i, word) in enumerate(firstline)
        index[i] = word
        data[word] = Vector{Float64}(undef, rows)
    end


    for row in 1:rows
        line = readline(io)
        try
            parse_split_line = parse.(Float64, split(line, '\t'))
            
            for (i, point) in enumerate(parse_split_line)
                data[index[i]][row] = point 
            end

        catch
            println("Couldn't parse line '", line, "'")
            continue
        end

    end
    
    close(io)
    return data
end

function get_febio_data_as_matrix(file::String)::Matrix{Float64}
    io = open(file, "r")
    
    rows = countlines(io) - 1; seek(io, 0)
    cols = length(split(readline(io))); seek(io, 0)

    data = Matrix{Float64}(undef, rows, cols)
    readline(io)

    for row in 1:rows
        data[row, :] = parse.(Float64, split(readline(io)))
    end
    

    close(io)
    return data
end


function get_node_positions(file::String)::Dict{String, Vector{Float64}}
    io = open(file, "r")

    node_pos = Dict{String, Vector{Float64}}()
    while eof(io) == false
        line = readline(io)
        
        if startswith(strip(line), "<Nodes name=")
            index = 1
            while true
                line = readline(io)
                r, node_id, x, y, z = @scanf(strip(line), "<node id=\"%i\">%f,%f,%f</node>", Int64, Float64, Float64, Float64)
                index += 1

                key = "N$node_id"
                node_pos[key] = [x, y, z]
                
                if r == 0
                    println("Succesfully read ", index, " lines")
                    break
                end
            end
        end
    end 
    close(io)
    return node_pos
end


function get_node_depths(sys::System, pos_dict::Dict{String, Vector{Float64}}, data_dict::Dict{String, Vector{Float64}})
    node_depths = Dict{String, Float64}()
    for (k, _) in data_dict
        if startswith(k, 'N') == false
            continue
        end
        node_depths[k] = sys.h - pos_dict[k][3]
    end
    return node_depths
end



function get_node_radius(pos_dict::Dict{String, Vector{Float64}}, data_dict::Dict{String, Vector{Float64}})
    node_depths = Dict{String, Float64}()
    for (k, _) in data_dict
        if startswith(k, 'N') == false
            continue
        end
        node_depths[k] = sqrt(pos_dict[k][1]^2 + pos_dict[k][2]^2)
    end
    return node_depths
end


