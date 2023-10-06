mutable struct SelfOrganizingMap
    data::Matrix{Float64}
    weights::Matrix{Float64}
    coords::Matrix{Float64}
    assignments::Vector{Int64}
    h₀::Function
    σ::Function
    function SelfOrganizingMap(data::AbstractMatrix{Float64}, n::Int64, h₀::R, σ::S) where {R<:Function, S<:Function}
        data = scaleData(data)
        weights = generateWeights(data, n)
        coords = generateCoords(n)
        assignments = Vector{Int64}(undef, size(data, 2))
        return new(data, weights, coords, assignments, h₀, σ)
    end
    function SelfOrganizingMap(data::Matrix{Float64}, weights::Matrix{Float64}, coords::Matrix{Float64}, assignments::Vector{Float64}, h₀::R, σ::S) where {R<:Function, S<:Function}
        return new(data, weights, coords, h₀, σ)
    end
end

function gauss(ri::Vector{Float64}, r0::Vector{Float64}, σ::Float64)
    return exp(-sqrt(sum((ri .- r0).^2)) / (2*σ^2))
end

function generateCoords(n)
    x = LinRange(-1.0, 1.0, n)' .* ones(n)
    y = ones(n)' .* LinRange(-1.0, 1.0, n)

    weights = Matrix{Float64}(undef, 2, n^2)
    Threads.@threads for i in 1:n^2
        weights[:,i] = [x[i], y[i]]
    end
    return weights
end

function generateWeights(D, n)
    return rand(-1.0:0.000001:1.0, size(D,1)-2, n)
end



function findBMU(SOM::SelfOrganizingMap, point::AbstractVector{Float64})
    distance = Vector{Float64}(undef, size(SOM.weights)[2])
    for (i, weight) in enumerate(eachcol(SOM.weights))
        distance[i] = sqrt(sum((weight .- point).^2))
    end
    return argmin(distance)
end

function updateWeights!(j::Int64, BMU_idx::Int64, SOM::SelfOrganizingMap, point::AbstractVector{Float64})
    for i in 1:size(SOM.weights)[2]
        SOM.weights[:,i] += SOM.h₀(j) .* gauss(SOM.coords[:,i], SOM.coords[:,BMU_idx], SOM.σ(j)) .* (point .- SOM.weights[:,i])
    end
end

function updateWeightsCluster!(j::Int64, BMU_idx::Int64, SOM::SelfOrganizingMap, point::AbstractVector{Float64})
    for i in 1:size(SOM.weights)[2]
        SOM.weights[:,i] += SOM.h₀(j) .* gauss(SOM.weights[:,i], SOM.weights[:,BMU_idx], SOM.σ(j)) .* (point .- SOM.weights[:,i])
    end
end


function shuffleColumns!(D)
    D = D[:,shuffle(1:size(D)[2])]
end

function train!(SOM, itrs::Int64)
    Threads.@threads for i in 1:itrs
        for point in eachcol(SOM.data[3:end,shuffle(1:size(SOM.data)[2])])
            BMU = findBMU(SOM, point)
            updateWeights!(i, BMU, SOM, point)
        end
    end
    getAssignments!(SOM)
end

function trainCluster!(SOM, itrs::Int64)
    Threads.@threads for i in 1:itrs
        for point in eachcol(SOM.data[:,shuffle(1:size(SOM.data)[2])])
            BMU = findBMU(SOM, point)
            updateWeightsCluster!(i, BMU, SOM, point)
        end
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



function plotData(SOM::SelfOrganizingMap, indices::AbstractVector, data::AbstractMatrix{Float64})
    out = assignData(SOM, data)
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


function getAssignments!(SOM::SelfOrganizingMap)
    for (i, point) in enumerate(eachcol(SOM.data[3:end,:]))
        SOM.assignments[i] = findBMU(SOM, point)
    end
end



function assignData(SOM::SelfOrganizingMap, data::AbstractMatrix)
    num = maximum(SOM.assignments)

    out = Vector{Matrix{Float64}}(undef, num)
    for i in 1:num
        out[i] = data[:,SOM.assignments.==i]
    end
    return out
end



function getCenterStatistics(SOM::SelfOrganizingMap, data::AbstractMatrix)
    set = assignData(SOM, data[3:end,:])

    mean_std = Vector{Tuple{Vector{Float64}, Vector{Float64}}}(undef, length(set))
    for (i, center) in enumerate(eachcol(SOM.weights))
        out = (vec(mean(set[i]; dims=2)), vec(std(set[i]; dims=2)))
        mean_std[i] = out
    end
    return mean_std
end


functio



function getBars(SOM::SelfOrganizingMap)

end






