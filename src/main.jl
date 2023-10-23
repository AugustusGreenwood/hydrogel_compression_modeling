


module DataInteract
    using Scanf, DataFrames, CSV, Statistics
    export getNodePositions, getElementConnectivity, getElementCenterOfMass, getData, getDataAtIndex, getAllDataAtIndex, scaleData, getTimepointData, getTimepointIndices
    include("FEBioInteract.jl")
    include("DataIO.jl")
end

module KMeans
    using Statistics, Plots, Clustering, DelimitedFiles, ..DataInteract
    export  splitDataByAssignment, plotKMeans, getCenterMeanStd, splitCenterStatsByVariable, 
            getAllResultsAtTimepoint, getAllCenterStatisticsAtTimepoint, getAllCenterStatisticsAtTimepoint,
            kmeans, outputAllCenterStatisticsAtTimepoint, plotSets, updateKMeans!
    include("KMeans.jl")
end

module SelfOrganizingMaps
    using Statistics, Plots, Random, ..DataInteract
    export getCenterMeanStd, splitDataByAssignment, plotSOM, train!, SelfOrganizingMap, splitCenterStatsByVariable, getAllMapsAtTimepoint
    include("SelfOrganizingMaps.jl")
end

module Waveform
    include("./CreateWaveform.jl")
end
