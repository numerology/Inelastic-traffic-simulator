function [ rate ] = maxmin(weights, wrapper_capacities, demands)
% Solve the maxmin rate allocation problem under the resource constraint.
% Cache is used to improve the simulation speed.
% --------------------------------
% Parameters:
% --------------------------------
% - userTypes : array of user types, 1 x nUsers matrix
% - weights: array of user weights, nUsers x 1 matrix
% - capacities: array of resources capacities, 1 x R matrix
% - demands: resource demands per user, R x nUsers matrix
% - cache: a map to store existing results
% --------------------------------
% Ret:
% --------------------------------
% - rate: rate allocation to each user, 1 x nUsers matrix
% - lambda: redundant, keep for backward comaptibility.
% First need to filter out the zeroes in the input
if max(weights) == 0
    rate = zeros(size(weights));
    return
end
nonZeroIdx = find(weights);

filteredWeights = weights(nonZeroIdx);
filteredDemands = demands(:, nonZeroIdx);    
filteredRate = maxminsolve(filteredWeights, ...
    wrapper_capacities, filteredDemands);
assert(all(filteredRate >= 0));
% change back
rate = zeros(size(weights))';
rate(nonZeroIdx) = filteredRate;
end