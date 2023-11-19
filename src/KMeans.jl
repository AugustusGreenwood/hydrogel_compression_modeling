

#=
This converts the vector in "getClusterPositionStatistics" so that either
the columns are cluster based (:column) or so that columns are variable
based (:variable)
=#
#=
For loop explanation:
input is a n length vector, where n is the number of clusters, of NxM matrices
where N = number of variables, M = 3 (mean, std, size of cluster)
We need to go from each vector being NxM to 1xN*M. reduce(hcat, eachrow(stat)')
accomplishes this.
Now we have a n length vector which is composed of 1xN*M matrics, so 
they just need to be stack from vector form into a complete matrix.
That what the final 'reduce(vcat, [...])' does. Don't @ me if it fails
=#
function collectClutserStatistics(stats::Vector{Matrix{Float64}}, type::Symbol)
    if type == :cluster
        return reduce(hcat, stats)
    elseif type == :variable
        return reduce(vcat, [reduce(hcat, eachrow(stat)') for stat in stats])
    else
        println("Didn't understand type, defaulting to cluster")
        return reduce(hcat, stats)
    end
end

#=
Gets the mean, std and N of the clusters belonging to a centroid
outputted this way for use in graphpad. should be then passed to
=#
function getClusterPositionStatistics(result::KmeansResult, data::Matrix{Float64})
    splits = splitDataByAssignments(result, data)

    stats = Vector{Matrix{Float64}}(undef, length(splits))
    for (i, split) in enumerate(splits)
        stats[i] = [mean(split, dims=2) std(split, dims=2) repeat([size(split, 2)], 7)]
    end
    return stats
end

#=
Gives a plot of a slice of the 5 dimentional data based on indices var.
plots the data and center together
=#
function plotKMeansResult(result::KmeansResult, data::Matrix{Float64}; indices::Vector{Int64}=[1, 2, 3])
    splits = splitDataByAssignments(result, data)
    colors = [:red, :blue, :green, :black, :purple]

    plot()
    for (i, split) in enumerate(splits)
        scatter!(split[indices[1], :], split[indices[2], :], split[indices[3], :], color=colors[i])
        scatter!(result.centers[indices[1], [i]], result.centers[indices[2], [i]], result.centers[indices[3], [i]], color=colors[i], markersize=15, markeralpha=0.5)
    end
    plot!()
end
#=
A kmeans result has integers on what cluster that data point is assigned to.
So we iterate through the number of clusters, get a bitwise vector to index
the data columns based on cluster number and return that
=#
function splitDataByAssignments(result::KmeansResult, data::Matrix{Float64})::Vector{Matrix{Float64}}
    clusters = size(result.centers, 2)

    out = Vector{Matrix{Float64}}(undef, clusters)
    for i in 1:clusters
        indices = result.assignments .== i
        out[i] = data[:, indices]
    end
    return out
end

#=
Functions that actually do the scaling. One replaces the data, very usefull if
I want to do something like a view since I don't usually want to change the
radius and height rows and I can pass the rest as a view, have it changed and
the data matrix is still like the original. The other returns a whole new matrix
=#
function scaleData!(data::AbstractMatrix{Float64}, fn::Function)
    for (i, row) in enumerate(eachrow(data))
        data[i, :] = fn(row)
    end
end
function scaleData(data::Matrix{Float64}, fn::Function)::Matrix{Float64}
    scaled_data = similar(data)
    for (i, row) in enumerate(eachrow(data))
        scaled_data[i, :] = fn(row)
    end
    return scaled_data
end

function _get_cluster_distance_statistics(centroid::Vector{Float64}, data::Matrix{Float64})
    distance = sqrt.(sum((centroid .- data) .^ 2, dims=1)) # Equclidian distance
    return mean(distance), std(distance)
end

#=
Scaling functions. Seems like scaling by the standard deviation is a good choice
as it perserves the sign, but we'll see
=#
function _0_to_1_scaling(row::AbstractVector)::Vector{Float64}
    min, max = minimum(row), maximum(row)
    return @. (row - min) / (max - min)
end

function _neg1_to_1_scaling(row::AbstractVector)::Vector{Float64}
    min, max = minimum(row), maximum(row)
    return @. 2 * (row - min) / (max - min) - 1
end

function _std_dev_scaling(row::AbstractVector)::Vector{Float64}
    return row ./ std(row)
end

function _0_to_1_scaling(ind_row::AbstractVector, whole_row::AbstractVector)::Vector{Float64}
    min, max = minimum(whole_row), maximum(whole_row)
    return @. (ind_row - min) / (max - min)
end

function _neg1_to_1_scaling(ind_row::AbstractVector, whole_row::AbstractVector)::Vector{Float64}
    min, max = minimum(whole_row), maximum(whole_row)
    return @. 2 * (ind_row - min) / (max - min) - 1
end

function _std_dev_scaling(ind_row::AbstractVector, whole_row::AbstractVector)::Vector{Float64}
    return ind_row ./ std(whole_row)
end
