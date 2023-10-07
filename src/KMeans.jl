function getCenterStatistics(result::KmeansResult, data::AbstractMatrix)
    set = assignData(result, data[3:end,:])

    mean_std = Vector{Tuple{Vector{Float64}, Vector{Float64}}}(undef, length(set))
    for (i, center) in enumerate(eachcol(result.centers))
        out = (vec(mean(set[i]; dims=2)), vec(std(set[i]; dims=2)))
        mean_std[i] = out
    end
    return mean_std
end









