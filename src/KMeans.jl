using Clustering

function splitDataByAssignment(result::KmeansResult, data::AbstractMatrix{Float64})
    num = maximum(result.assignments)

    out = Vector{Matrix{Float64}}(undef, num)
    for i in 1:num
        out[i] = data[:,result.assignments.==i]
    end
    return out
end

function plotKMeans(result::KmeansResult, indices::AbstractVector, data::AbstractMatrix{Float64})
    out = splitDataByAssignment(result, data[3:end, :])
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


function getCenterMeanStd(result::KmeansResult, data::AbstractMatrix)
    set = splitDataByAssignment(result, data[3:end,:])

    mean_std = Vector{Matrix{Float64}}(undef, length(set))
    for i in eachindex(set)
        mean_std[i] = [vec(mean(set[i]; dims=2)) vec(std(set[i]; dims=2))]
    end
    return reduce(hcat, mean_std)
end

function splitCenterStatsByCluster(stats)
    by_cluster = reduce(vcat, [reduce(hcat, [stats[i,2*j-1:2*j] for i in 1:5]') for j in 1:4])
    return by_cluster
end

function getAllResultsAtTimepoint(dir::String, k::Int64)
    data = getTimepointData(dir)
    scaled_data = [scaleData(d) for d in data]
    results = [(kmeans(d[3:end,:], k), d) for d in scaled_data]

    return results
end

function getAllCenterStatisticsAtTimepoint(dir::String, k::Int64; sample::Bool=false, indices::AbstractVector=[1,2,3])
    results = getAllResultsAtTimepoint(dir, k)
    stats = [getCenterMeanStd(d[1], d[2]) for d in results]

    if sample == true
        for result in results
            plotKMeans(result[1], indices, result[2])
        end
    end
    return stats
end


function outputAllCenterStatisticsAtTimepoint(dir::String, k::Int64; by_cluster=true, sample::Bool=false, indices::AbstractVector=[2,3,4])
    output_names_cluster = joinpath.(dir, ["eor_kmeans_cluster.csv", "eor_dwell_kmeans_cluster.csv", "eoc_kmeans_cluster.csv", "eoc_dwell_kmeans_cluster.csv"])
    output_names_variable = joinpath.(dir, ["eor_kmeans_var.csv", "eor_dwell_kmeans_var.csv", "eoc_kmeans_var.csv", "eoc_dwell_kmeans_var.csv"])

    stats = getAllCenterStatisticsAtTimepoint(dir, k; sample=sample, indices=indices)

    if by_cluster == true
        stats = splitCenterStatsByCluster.(stats)
        output_names = output_names_cluster
    elseif by_cluster == false
        output_names = output_names_variable
    elseif by_cluster == :both
        cluster_stats = splitCenterStatsByCluster.(stats)
        var_stats = stats
        for i in eachindex(stats)
            writedlm(output_names_variable[i], var_stats[i], ',')
            writedlm(output_names_cluster[i], cluster_stats[i], ',')
        end
        return
    end

    for (i, stat) in enumerate(stats)
        writedlm(output_names[i], stat, ',')
    end
end

function getSets(dir::String, k::Int64, tol::Float64, timepoint_index::Int64)
    # Open to different options, but right now the sets are plotted by the radial center of mass
    sortFunc(x) = sum(x[1,:]) / size(x, 2) 

    data = getTimepointData(dir)[timepoint_index]
    scaled_data = scaleData(data)
    result = kmeans(scaled_data[3:end,:], k; tol=tol)
    sets = sort(splitDataByAssignment(result, data[1:2,:]), by=sortFunc)
    return sets
end



function plotSets(dir::String, k::Int64, tol::Float64, timepoint_index::Int64; colors=[:red, :blue, :green, :purple, :black, :yellow], plot_options...)
    # Open to different options, but right now the sets are plotted by the radial center of mass
    sortFunc(x) = sum(x[1,:]) / size(x, 2) 

    data = getTimepointData(dir)[timepoint_index]
    scaled_data = scaleData(data)
    result = kmeans(scaled_data[3:end,:], k; tol=tol)
    sets = sort(splitDataByAssignment(result, data[1:2,:]), by=sortFunc)

    plts = plot()
    for (i, set) in enumerate(sets)
        scatter!(set[1,:], set[2,:], color=colors[i]; plot_options...);
    end
    return plts
end

function plotSets(result::KmeansResult, data::AbstractMatrix; colors=[:red, :blue, :green, :purple, :black, :yellow], plot_options...)
    # Open to different options, but right now the sets are plotted by the radial center of mass
    sortFunc(x) = sum(x[1,:]) / size(x, 2) 

    sets = sort(splitDataByAssignment(result, data[1:2,:]), by=sortFunc)

    plts = plot()
    for (i, set) in enumerate(sets)
        scatter!(set[1,:], set[2,:], color=colors[i]; plot_options...);
    end
    return plts
end

function plotSets!(plt::AbstractPlot, dir::String, k::Int64, tol::Float64, timepoint_index::Int64; colors=[:red, :blue, :green, :purple, :black, :yellow], plot_options...)
    # Open to different options, but right now the sets are plotted by the radial center of mass
    sortFunc(x) = sum(x[1,:]) / size(x, 2) 

    data = getTimepointData(dir)[timepoint_index]
    scaled_data = scaleData(data)
    result = kmeans(scaled_data[3:end,:], k; tol=tol)
    sets = sort(splitDataByAssignment(result, data[1:2,:]), by=sortFunc)

    for (i, set) in enumerate(sets)
        scatter!(set[1,:], set[2,:], color=colors[i]; plot_options...);
    end
end

function updateKMeans!(result::KmeansResult, data::AbstractMatrix, k::Int64; kmeans_options...)
    bmus = __getAllBestMatchingUnits(result, data)
    result = kmeans(data, k, init=bmus; kmeans_options...)
end

function __getBestMatchingUnitIndex(data::AbstractMatrix, center::AbstractVector{Float64})
    distance = Vector{Float64}(undef, size(data, 2))
    for (i, point) in enumerate(eachcol(data))
        distance[i] = sqrt(sum((center .- point) .^ 2))
    end
    return argmin(distance)
end


function __getAllBestMatchingUnits(result::KmeansResult, data::AbstractMatrix)
    bmus = Vector{Int64}(undef, size(result.centers, 2))
    for (i, center) in enumerate(eachcol(result.centers))
        bmus[i] = __getBestMatchingUnitIndex(data, center)
    end
    return bmus
end
