


module DataInteract
    using Scanf, DataFrames, CSV, Statistics
    export getNodePositions, getElementConnectivity, getElementCenterOfMass, getData, getDataAtIndex, getAllDataAtIndex, scaleData
    include("FEBioInteract.jl")
    include("DataIO.jl")
end

module KMeans
    using Statistics, Plots, ..DataInteract
    plotlyjs()
    export getCenterStatistics, assignData, plotKMeans
    include("KMeans.jl")
end

module SelfOrganizingMaps
    using Statistics, Plots, Random, ..DataInteract
    plotlyjs()
    export getCenterMeanStd, splitDataByAssignment, plotSOM, train!, SelfOrganizingMap
    include("SOM.jl")
end

module Waveform
    include("./CreateWaveform.jl")
end
