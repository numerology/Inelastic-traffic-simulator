function [inverseP, inverseLevel] = inversewaterfill(total, ceiling)
% inversewaterfill Inverse water filling function. Try to solve the minmax
% problem of the amount of water in each bucket, subject to the constraint that
% the height of each bucket is upperbounded by vector ceiling.

L = length(ceiling); % number of buckets
totalVol = sum(ceiling);
assert(total <= totalVol, 'No feasible allocation found for water fill.');

[p, level] = waterfill(max(ceiling) - ceiling, totalVol - total);
inverseLevel = max(ceiling) - level;
inverseP = ceiling - p;
end