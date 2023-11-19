include("ClusterAnalysis.jl")
using Statistics
import .ClusterAnalysis: readFEBioDataFile, getAllFEBioDataAtTime, collectFEBioData,
    getFEBioDataAtLocation, getFEBioDataAtTime, writeClusterStatistics,
    writeClusterStatistics, writeClusterStatisticsForProject,
    writeClusterStatisticsForDirectory, collectClutserStatistics,
    getClusterPositionStatistics, plotKMeansResult,
    splitDataByAssignments, scaleData!, scaleData, kmeans, kmeans!,
    _05_times, _075_times, _1_times, _sin_times, readAllProjectData,
    _0_to_1_scaling, _std_dev_scaling, _neg1_to_1_scaling

mat_data_1, mat_data_075, mat_data_05, mat_data_sin = readAllProjectData(_1_times[:EOR], _075_times[:EOR], _05_times[:EOR], _sin_times[:EOR])
data = [view(mat_data_1, 3:7, :), view(mat_data_075, 3:7,:), 
        view(mat_data_05, 3:7,:), view(mat_data_sin, 3:7,:)]

function scale!(a::Vector, fn::Function)
    display(typeof(a))
    for (i, row) in enumerate(eachrow(hcat(a...)))
        for j in eachindex(a)
            a[j][i,:] = fn(a[j][i,:], row)
        end
    end
end

function _0_to_1_scaling(ind_row::AbstractVector, whole_row::AbstractVector)::Vector{Float64}
    min, max = minimum(whole_row), maximum(whole_row)
    return @. (ind_row - min) / (max - min)
end

function _neg1_to_1_scaling(ind_row::AbstractVector, whole_row::AbstractVector)::Vector{Float64}
    min, max = minimum(whole_row), maximum(whole_row)
    return @. 2 * (ind_row - min) / (max - min) - 1
end

function _std_dev_scaling(ind_row::AbstractVector, whole_row::AbstractVector)::Vector{Float64}
    return ind_row ./ std(whole_row)
end

function writeClusterStatistics(data::Matrix{Float64}, out_file_name::String, var_clust::Symbol)
    result = kmeans(view(data, 3:7, :), 3; display=:final, tol=0.0)
    vec_stats = getClusterPositionStatistics(result, data)
    mat_stats = collectClutserStatistics(vec_stats, var_clust)
    writeClusterStatistics(mat_stats, out_file_name)
end


time = :EOR;
time_string = "EOR";
scaling = _neg1_to_1_scaling;
scale_string = "n11";
mat_data_1 = collectFEBioData(getAllFEBioDataAtTime("1s-dwell_3s-period/", _1_times[time]), true)
mat_data_075 = collectFEBioData(getAllFEBioDataAtTime("0.75s-dwell_3s-period/", _075_times[time]), true);
mat_data_05 = collectFEBioData(getAllFEBioDataAtTime("0.5s-dwell_3s-period/", _05_times[time]), true);
mat_data_sin = collectFEBioData(getAllFEBioDataAtTime("sinusoid/", _sin_times[time]), true);
data = [view(mat_data_1, 3:7, :), view(mat_data_075, 3:7,:), 
        view(mat_data_05, 3:7,:)];
scale!(data, scaling);

writeClusterStatistics(mat_data_1,         "1s-dwell_3s-period/scaled_across/cluster/$(time_string)_$(scale_string)_cluster_whole_scale.csv", :cluster);
writeClusterStatistics(mat_data_075,    "0.75s-dwell_3s-period/scaled_across/cluster/$(time_string)_$(scale_string)_cluster_whole_scale.csv", :cluster);
writeClusterStatistics(mat_data_05,      "0.5s-dwell_3s-period/scaled_across/cluster/$(time_string)_$(scale_string)_cluster_whole_scale.csv", :cluster);
writeClusterStatistics(mat_data_sin,                 "sinusoid/scaled_across/cluster/$(time_string)_$(scale_string)_cluster_whole_scale.csv", :cluster);

writeClusterStatistics(mat_data_1,         "1s-dwell_3s-period/scaled_across/variable/$(time_string)_$(scale_string)_variable_whole_scale.csv", :variable);
writeClusterStatistics(mat_data_075,    "0.75s-dwell_3s-period/scaled_across/variable/$(time_string)_$(scale_string)_variable_whole_scale.csv", :variable);
writeClusterStatistics(mat_data_05,      "0.5s-dwell_3s-period/scaled_across/variable/$(time_string)_$(scale_string)_variable_whole_scale.csv", :variable);
writeClusterStatistics(mat_data_sin,                 "sinusoid/scaled_across/variable/$(time_string)_$(scale_string)_variable_whole_scale.csv", :variable);
display(mat_data_1)


