global const EOC_INDEX::Int64, EOC_TIME::Float64, EOR_INDEX::Int64, EOR_TIME::Float64 = __get_eoc_eor_locations()::Tuple{Int64, Float64, Int64, Float64}

function __get_eoc_eor_locations()::Tuple{Int64, Float64, Int64, Float64}
    file = open("locations.dat", "r")
    eoc_line, eor_line = readline(file), readline(file)

    r, eoc_indx, eoc_time = @scanf(eoc_line, "%i,%f", Int64, Float64)
    r, eor_indx, eor_time = @scanf(eor_line, "%i,%f", Int64, Float64)

    return eoc_indx, eoc_time, eor_indx, eor_time
end

function __get_r_h(feb_file::String)
    node_loc = getFEBioNodePositions(feb_file)

    node_rh = Dict{Int64, Matrix{Float64}}()
    for (node, loc) in node_loc
        r = sqrt(loc[1]^2 + loc[2]^2)
        h = loc[3]
        node_rh[node] = [r h]
    end
    return node_rh
end


function __get_r_h_data(data_file::String, index::Int64)
    data = getFEBioData(data_file; out_type=:dict)
    
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
    data = getFEBioData(data_file; out_type=:dict)
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

function getAllData(root::String, index::Int64)
    data_files = data_files = [file for file in readdir("$root/", join=true) if endswith(file, ".dat")]
    data_files = data_files[[7,2,3,4,5,6,1]]
    vec_data = getFEBioDataAtIndex.(data_files, index)
    data = reduce(hcat, vec_data)'
    return data
end

function outputVectorData(data_file::String, feb_file::String, out_file::String; seperate::Bool=false)
	rhd_dict = __get_r_h_eoc_eor_data(data_file, feb_file)
	rhd_mat = __dict_to_mat_sort(rhd_dict)
    if seperate == false
        writedlm(out_file, rhd_mat', ',')
    else
        r_eoc_mat = [rhd_mat[:,1] rhd_mat[:,3]]
        h_eoc_mat = [rhd_mat[:,2] rhd_mat[:,3]]
        r_eor_mat = [rhd_mat[:,1] rhd_mat[:,4]]
        h_eor_mat = [rhd_mat[:,2] rhd_mat[:,4]]

        writedlm("r_eoc_$out_file", r_eoc_mat, ',')
        writedlm("h_eoc_$out_file", h_eoc_mat, ',')
        writedlm("r_eor_$out_file", r_eor_mat, ',')
        writedlm("h_eor_$out_file", h_eor_mat, ',')
    end
end



function outputMatrixData(data_file::String, feb_file::String, out_file::String)
	rhd_dict = __get_r_h_eoc_eor_data(data_file, feb_file)
	rhd_mat = __dict_to_mat_sort(rhd_dict)
    r = unique(rhd_mat[:,1])
    h = unique(rhd_mat[:,2])
    eoc = reshape(rhd_mat[:,3], 26, 26)
    eor = reshape(rhd_mat[:,4], 26, 26)
    eoc_mat = hcat([0;h], vcat(r', eoc))
    eor_mat = hcat([0;h], vcat(r', eor))
    writedlm("eoc_$out_file", eoc_mat, ',')
    writedlm("eor_$out_file", eor_mat, ',')
end





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


function assignData(result::KmeansResult, data::AbstractMatrix{Float64})
    num = maximum(result.assignments)

    out = Vector{Matrix{Float64}}(undef, num)
    for i in 1:num
        out[i] = data[:,result.assignments.==i]
    end
    return out
end


function assignData(result::KmedoidsResult, data::AbstractMatrix{Float64})
    num = maximum(result.assignments)

    out = Vector{Matrix{Float64}}(undef, num)
    for i in 1:num
        out[i] = data[:,result.assignments.==i]
    end
    return out
end


function plotData(result::KmeansResult, data::AbstractMatrix{Float64})
    out = assignData(result, data)

    if size(data, 1) == 3 
        plt = scatter(size=(900,700))
        for (i, set) in enumerate(out)
            scatter!(set[1,:], set[2,:], set[3,:], label="cluster $i")
            scatter!([result.centers[1,i]], [result.centers[2,i]], [result.centers[3,i]], color=plt.series_list[end][:linecolor], label="center $i", markersize=20, markeralpha=0.5)
        end
        display(plt)
    else
        plt = scatter(size=(900,700))
        for (i, set) in enumerate(out)
            scatter!(set[1,:], set[2,:], label="cluster $i")
            scatter!([result.centers[1,i]], [result.centers[2,i]], color=plt.series_list[end][:linecolor], label="center $i", markersize=20, markeralpha=0.5)
        end
        display(plt)
    end
    return plt
end



function plotData(result::KmedoidsResult, index::Int64, data::AbstractMatrix{Float64})
    out = assignData(result, data)
    labels = ["rpos", "hpos", "pressure", "strain 1st invariant", "strain 2nd invariant", "strain 3rd invariant", "fluid flux magnitude"]

    plt = scatter(size=(900,700))
    for (i, set) in enumerate(out)
        to_plot = [set[1:2,:]; set[index,:]']
        scatter!(to_plot[1,:], to_plot[2,:], to_plot[3,:], label="cluster $i")
        scatter!([result.centers[1,i]], [result.centers[2,i]], [result.centers[index,i]], color=plt.series_list[end][:linecolor], label="center $i", markersize=20, markeralpha=0.5)
        zlabel!(labels[index])
        ylabel!(labels[2])
        xlabel!(labels[1])
    end
    display(plt)
    return plt
end

function plotData(result::KmeansResult, indices::AbstractVector, data::AbstractMatrix{Float64})
    out = assignData(result, data[3:end, :])
    labels = ["pressure", "strain 1st invariant", "strain 2nd invariant", "strain 3rd invariant", "fluid flux magnitude"]

    plt = scatter(size=(900,700))
    for (i, set) in enumerate(out)
        data_to_plot = set[indices, :]
        center_to_plot =  result.centers[indices, i]
        scatter!(data_to_plot[1,:], data_to_plot[2,:], data_to_plot[3,:])
        scatter!(center_to_plot[1,:], center_to_plot[2,:], center_to_plot[3,:], markersize=15, markeralpha=0.5, color=plt.series_list[end][:linecolor]) 
        xlabel!(labels[indices[1]])
        ylabel!(labels[indices[2]])
        zlabel!(labels[indices[3]])
    end
    display(plt)
    return plt
end

function getDistanceMatrix(mat)
    n = size(mat)[2]
    distance_matrix = zeros(Float64, n, n) 
    for (i, wi) in enumerate(eachcol(mat))
        for (j, wj) in enumerate(eachcol(mat))
            if j > i
                distance = sqrt(sum((wi .- wj).^2))
                distance_matrix[j,i] = distance_matrix[i,j] = distance
            end
        end
    end
    
    return distance_matrix / maximum(distance_matrix)
end


function getSimilarMatrix(mat::AbstractMatrix, dim::Int64)
    n,d = size(mat)

    out = Matrix{Float64}(undef, n, dim)
    for (i, row) in enumerate(eachrow(mat))
        out[i, :] = rand(minimum(row):1e-15:maximum(row), dim)
    end
    return out
end


function getHopkinsStatistic(X::AbstractMatrix{Float64}, m::Int64)
    d,n = size(X)
    X_tilde = X[:, randperm(n)[1:m]]
    Y = getSimilarMatrix(X, m)
    temp = Vector{Float64}(undef, n)

    ui = Vector{Float64}(undef, m)
    for (i, Yi) in enumerate(eachcol(Y))
        for (j, Xi) in enumerate(eachcol(X))
            temp[j] = sqrt(sum( (Yi .- Xi).^2))
        end
        ui[i] = minimum(temp)
    end

    wi = Vector{Float64}(undef, m)
    for (i, Xi) in enumerate(eachcol(X_tilde))
        for (j, Xj) in enumerate(eachcol(X))
            if Xi != Xj
                temp[j] = sqrt(sum( (Xi .- Xj).^2 ))
            end
        end
        wi[i] = minimum(temp)
    end

   u, w = sum(ui), sum(wi)

    return u / (u + w)
end


function getAverageHopkinsStatistic(X::AbstractMatrix{Float64})
    d,n = size(X)

    H_avg = Vector{Float64}(undef, n - 2)
    Threads.@threads for c in 2:n-1
        H_avg[c-1] = getHopkinsStatistic(X, c)
    end
    return sum(H_avg) / length(H_avg)
end

function getAllData(path::String, index::Int64)
    data_files = [file for file in readdir(path, join=true) if endswith(file, ".dat")]
    data_files = data_files[[7,2,3,4,5,6,1]]
    vec_data = getFEBioDataAtIndex.(data_files, index)
    return reduce(hcat, vec_data)'
end




function getCenterStatistics(result::KmeansResult, data::AbstractMatrix)
    set = assignData(result, data[3:end,:])

    mean_std = Vector{Tuple{Vector{Float64}, Vector{Float64}}}(undef, length(set))
    for (i, center) in enumerate(eachcol(result.centers))
        out = (vec(mean(set[i]; dims=2)), vec(std(set[i]; dims=2)))
        mean_std[i] = out
    end
    return mean_std
end


