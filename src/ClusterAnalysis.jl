module ClusterAnalysis

include("DataIO.jl")
export readFEBioDataFile, getAllFEBioDataAtTime, collectFEBioData
export getFEBioDataAtLocation, getFEBioDataAtTime, writeClusterStatistics
export writeClusterStatistics, writeClusterStatisticsForProject
export writeClusterStatisticsForDirectory





using Clustering, Statistics, Plots
include("KMeans.jl")
export collectClutserStatistics, getClusterPositionStatistics, plotKMeansResult
export splitDataByAssignments, scaleData!, scaleData, kmeans, kmeans!
export _0_to_1_scaling, _std_dev_scaling, _neg1_to_1_scaling

# Something would be better (probably) some other way, but globals
# are the easiest by far

# end of relaxation,    
# end of relaxation dwell, 
# end of compression,      
# end of compression dwell
const global _05_times::Dict{Symbol,Float64} = Dict{Symbol,Float64}(
    :EOR => 601.0, :EORD => 601.5, :EOC => 602.5, :EOCD => 603.0
)
const global _075_times::Dict{Symbol,Float64} = Dict{Symbol,Float64}(
    :EOR => 594.897, :EORD => 596.25, :EOC => 596.25, :EOCD => 597.0
)
const global _1_times::Dict{Symbol,Float64} = Dict{Symbol,Float64}(
    :EOR => 594.546, :EORD => 595.446, :EOC => 596.046, :EOCD => 596.946
)
const global _sin_times::Dict{Symbol,Float64} = Dict{Symbol,Float64}(
    :EOR => 595.437, :EOC => 596.937
)

export _sin_times, _05_times, _075_times, _1_times



end


