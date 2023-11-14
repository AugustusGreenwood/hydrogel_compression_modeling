module ClusterAnalysis

include("DataIO.jl")
export readFEBioDataFile, getAllFEBioDataAtTime, collectFEBioData,
    getFEBioDataAtLocation, getFEBioDataAtTime


# Something would be better (probably) some other way, but globals
# are the easiest by far

# end of relaxation,    
# end of relaxation dwell, 
# end of compression,      
# end of compression dwell
global const _05_times::Dict{Symbol, Float64} = Dict{Symbol, Float64}(
    :EOR => 601.0, :EORD => 601.5, :EOC => 602.5, :EOCD => 603.0
)
global const _075_times::Dict{Symbol, Float64} = Dict{Symbol, Float64}(
    :EOR => 594.897, :EORD => 596.25, :EOC => 596.25, :EOCD => 597.0
)
global const _1_times::Dict{Symbol, Float64} = Dict{Symbol, Float64}(
    :EOR => 594.546, :EORD => 595.446, :EOC => 596.046, :EOCD => 596.946
)
global const _sin_times::Dict{Symbol, Float64} = Dict{Symbol, Float64}(
    :EOR => 595.437, :EOC => 596.937
)


export _sin_times, _05_times, _075_times, _1_times

using Clustering, Statistics, Plots
include("KMeans.jl")
export collectClutserStatistics, getClusterPositionStatistics, plotKMeansResult,
        splitDataByAssignments, scaleData!, scaleData, kmeans, kmeans!
        _0_to_1_scaling, _std_dev_scaling, _neg1_to_1_scaling
end
