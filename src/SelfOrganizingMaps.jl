mutable struct SelfOrganizingMap
    data::Matrix{Float64}
    weights::Matrix{Float64}
    coords::Matrix{Float64}
    assignments::Vector{Int64}
    h₀::Function
    σ::Function
    function SelfOrganizingMap(data::AbstractMatrix{Float64}, n::Int64, h₀::R, σ::S) where {R<:Function, S<:Function}
        data = scaleData(data)
        weights = __generateWeights(data, n)
        coords = __generateCoordinates(n)
        assignments = Vector{Int64}(undef, size(data, 2))
        return new(data, weights, coords, assignments, h₀, σ)
    end
end


function train!(SOM, itrs::Int64)
    Threads.@threads for i in 1:itrs
        for point in eachcol(__shuffleColumns(SOM.data[3:end,:]))
            BMU = __getBestMatchingUnitIndex(SOM, point)
            __updateWeights!(i, BMU, SOM, point)
        end
    end
    __getAssignments!(SOM)
end

function splitDataByAssignment(SOM::SelfOrganizingMap, data::AbstractMatrix)
    num = maximum(SOM.assignments)

    out = Vector{Matrix{Float64}}(undef, num)
    for i in 1:num
        out[i] = data[:,SOM.assignments.==i]
    end
    return out
end

function plotSOM(SOM::SelfOrganizingMap, indices::AbstractVector, data::AbstractMatrix{Float64})
    out = splitDataByAssignment(SOM, data)
    labels = ["pressure", "strain 1st invariant", "strain 2nd invariant", "strain 3rd invariant", "fluid flux magnitude"]

    plt = scatter(size=(900,700))
    for (i, set) in enumerate(out)
        data_to_plot = set[indices, :]
        center_to_plot =  SOM.weights[indices, i]
        scatter!(data_to_plot[1,:], data_to_plot[2,:], data_to_plot[3,:])
        scatter!(center_to_plot[1,:], center_to_plot[2,:], center_to_plot[3,:], markersize=15, markeralpha=0.5, color=plt.series_list[end][:linecolor]) 
        xlabel!(labels[indices[1]])
        ylabel!(labels[indices[2]])
        zlabel!(labels[indices[3]])
    end
    display(plt)
    return plt
end

function getCenterMeanStd(SOM::SelfOrganizingMap, data::AbstractMatrix)
    set = splitDataByAssignment(SOM, data[3:end,:])

    mean_std = Vector{Matrix{Float64}}(undef, length(set))
    for (i, center) in enumerate(eachcol(SOM.weights))
        mean_std[i] = [vec(mean(set[i]; dims=2)) vec(std(set[i]; dims=2))]
    end
    return reduce(hcat, mean_std)
end


function splitCenterStatsByVariable(SOM::SelfOrganizingMap, data::AbstractMatrix)
    stats = getCenterMeanStd(SOM, data)
    var_stats = [ reduce(hcat, [[stats[i][1][j], stats[i][2][j]] for i in 1:4]) for j in 1:5]
    return var_stats
end


function __gaussian(ri::Vector{Float64}, r0::Vector{Float64}, σ::Float64)
    return exp(-sqrt(sum((ri .- r0).^2)) / (2*σ^2))
end

function __generateCoordinates(n)
    x = LinRange(-1.0, 1.0, n)' .* ones(n)
    y = ones(n)' .* LinRange(-1.0, 1.0, n)

    weights = Matrix{Float64}(undef, 2, n^2)
    Threads.@threads for i in 1:n^2
        weights[:,i] = [x[i], y[i]]
    end
    return weights
end

function __generateWeights(D, n)
    return rand(-1.0:0.000001:1.0, size(D,1)-2, n)
end

function __getBestMatchingUnitIndex(SOM::SelfOrganizingMap, point::AbstractVector{Float64})
    distance = Vector{Float64}(undef, size(SOM.weights)[2])
    for (i, weight) in enumerate(eachcol(SOM.weights))
        distance[i] = sqrt(sum((weight .- point).^2))
    end
    return argmin(distance)
end

function __updateWeights!(j::Int64, BMU_idx::Int64, SOM::SelfOrganizingMap, point::AbstractVector{Float64})
    for i in 1:size(SOM.weights)[2]
        SOM.weights[:,i] += SOM.h₀(j) .* __gaussian(SOM.coords[:,i], SOM.coords[:,BMU_idx], SOM.σ(j)) .* (point .- SOM.weights[:,i])
    end
end

function __shuffleColumns(D)
    return D[:,shuffle(1:size(D)[2])]
end


function __getAssignments!(SOM::SelfOrganizingMap)
    for (i, point) in enumerate(eachcol(SOM.data[3:end,:]))
        SOM.assignments[i] = __getBestMatchingUnitIndex(SOM, point)
    end
end


