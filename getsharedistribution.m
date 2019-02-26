function [shareDist] = getsharedistribution(opSettings, loadDist)
% Return the share distribution aligned with each slice's load distribution
% Currently, only allows equal share vector.
% Params:
%   opSettings:
%   loadDist: V x B, relative load distribution of each slice
% Returns:
%   shareDist: V x B, share distribution across BSs of each slice
shareVec = opSettings.s_o;
V = size(loadDist, 1);
nBS = size(loadDist, 2);
% First need to renormalize per BS
normalizedLoadDist = loadDist;
for b = 1:nBS
    normalizedLoadDist(:, b) = (shareVec .* normalizedLoadDist(:, b)) ...
        / sum(shareVec .* normalizedLoadDist(:, b));
end
shareDist = normalizedLoadDist;
%{
for v = 1:V
    shareDist(v, :) = shareVec(v) * (shareDist(v, :)/nansum(shareDist(v, :)));
end
%}
shareDist(isnan(shareDist)) = 0;
end

