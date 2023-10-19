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


function splitCenterStatsByVariable(result::KmeansResult, data::AbstractMatrix)
    stats = getCenterMeanStd(result, data)
    var_stats = [ reduce(hcat, [[stats[i][1][j], stats[i][2][j]] for i in 1:4]) for j in 1:5]
    return var_stats
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


function outputAllCenterStatisticsAtTimepoint(dir::String, k::Int64; sample::Bool=false, indices::AbstractVector=[2,3,4])
    var_labels = ["pressure", "strain 1st invariant", "strain 2nd invariant", "strain 3rd invariant", "fluid flux magnitude"]
    cluster_labels = [0, "cluster 1 mean", "cluster 1 mean", "cluster 2 mean", "cluster 2 std", "cluster 3 mean", "cluster 3 std", "cluster 4 mean", "cluster 4 std"]
    output_names = joinpath.(dir, ["eor_kmeans.csv", "eor_dwell_kmeans.csv", "eoc_kmeans.csv", "eoc_dwell_kmeans.csv"])

    to_output = Matrix{Any}(undef, 9, 6)
    to_output[:,1] = cluster_labels
    to_output[1,2:end] = var_labels

    stats = getAllCenterStatisticsAtTimepoint(dir, k; sample=sample, indices=indices)
    for (i, stat) in enumerate(stats)
        to_output[2:end,2:end] = stat'
        writedlm(output_names[i], to_output, ',')
    end
end


function plotSets(dir::String, k::Int64, timepoint_index::Int64, options...)
    data = getTimepointData(dir)[timepoint_index]
    scaled_data = scaleData(data)
    result = kmeans(scaled_data[3:end,:], 4)
    sets = splitDataByAssignment(result, data[1:2,:])

    plots = plot(legend=false, size=(800,500), xlims=(0,3), ylims=(0,1.6), options...)
    for set in sets
        scatter!(set[1,:], set[2,:])
    end
    plot!()
    return plots
end



