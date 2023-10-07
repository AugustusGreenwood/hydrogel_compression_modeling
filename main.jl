using Plots; plotlyjs()
using Random
using Scanf
using DataFrames
using DelimitedFiles
using Statistics
using CSV
using Clustering



module FEBioInteract
export getNodePositions, getElementConnectivity, getElementCenterOfMass, getData, getDataAtIndex, getAllDataAtIndex
    include("./FEBioInteract.jl")
end

module DataAnalysis
    include("./DataAnalysis.jl")
end

module Waveform
    include("./CreateWaveform.jl")
end

module Clustering
    include("./Clustering.jl")
end

