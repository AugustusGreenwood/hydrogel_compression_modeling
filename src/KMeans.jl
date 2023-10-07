using Clustering

function getCenterStatistics(result::KmeansResult, data::AbstractMatrix)
    set = assignData(result, data[3:end,:])

    mean_std = Vector{Tuple{Vector{Float64}, Vector{Float64}}}(undef, length(set))
    for (i, center) in enumerate(eachcol(result.centers))
        out = (vec(mean(set[i]; dims=2)), vec(std(set[i]; dims=2)))
        mean_std[i] = out
    end
    return mean_std
end

function assignData(result::KmeansResult, data::AbstractMatrix{Float64})
    num = maximum(result.assignments)

    out = Vector{Matrix{Float64}}(undef, num)
    for i in 1:num
        out[i] = data[:,result.assignments.==i]
    end
    return out
end

function plotKMeans(result::KmeansResult, indices::AbstractVector, data::AbstractMatrix{Float64})
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







