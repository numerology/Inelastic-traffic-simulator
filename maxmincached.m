function [ rate, lambda ] = maxmincached(userTypes, weights, capacities, ...
    demands, cache)
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
nonZeroIdx = find(weights);
filteredWeights = weights(nonZeroIdx);
% Adjust the order according to the order of user_weights
[orderedFilteredWeights, weightsIndex] = sort(filteredWeights);
[tmp, sortBack] = sort(weightsIndex);
sortedDemands = demands(:, weightsIndex);

T = max(userTypes);
nUserTypes = zeros(1,T);
for t = userTypes
    nUserTypes(t) = nUserTypes(t) + 1;
end

keyGenerated = keygen(nUserTypes);
if(cache.isKey(keyGenerated))
    orderedFilteredRate = cache(keyGenerated);
    lambda = -1;
    %disp('cache hit');
else
    orderedFilteredRate = maxminsolve(orderedFilteredWeights, ...
        capacities, sortedDemands);
    
    assert(all(orderedFilteredRate >= 0));
end
% change back to the original order
filtered_rate = orderedFilteredRate(sortBack);
rate = zeros(size(weights))';
rate(nonZeroIdx) = filtered_rate;

lambda = 1; % REDUNDANT, but want to mimic the behavior of fminimax for 
% compatibility.
end

