include("ClusterAnalysis.jl")
import .ClusterAnalysis: readFEBioDataFile, getAllFEBioDataAtTime, 
    collectFEBioData, getFEBioDataAtLocation, getFEBioDataAtTime,
    collectClutserStatistics, getClusterPositionStatistics, plotKMeansResult,
    splitDataByAssignments, scaleData!, scaleData,
    _0_to_1_scaling, _std_dev_scaling, _neg1_to_1_scaling;




