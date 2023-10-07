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
    vec_data = getDataAtIndex.(data_files, index)
    return reduce(hcat, vec_data)'
end


