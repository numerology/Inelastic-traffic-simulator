function [loadDist, mask] = getloaddistribution(opSettings, netSettings, bsAssociation, ...
    simulationTime)
% Get the relative load distribution across base stations per slice.
% Params:
%   opSettings: set up
%   bsAssociation: bsAssociation(u, t) = the index of the BS user u is
%   associated with at time t
%   simulationTime: duration of the simulation
% Returns:
%   loadDist: V x bsNs, loadDist(v, b) = mean number of users associated with
%   base station b through the simulation
%   mask: 1 x bsNs, 1 if the BS has a nonzero load, 0 if the BS has 0 load
%   thus can be excluded from the discussion.

V = opSettings.operators;
nBS = netSettings.bsNS;
opBelongs = opSettings.ops_belongs;
loadDist = zeros(V, nBS);
mask = zeros(1, nBS);
for t = 1:simulationTime
    for v = 1:V
        for b = 1:nBS
            loadDist(v, b) = loadDist(v, b) ...
                + sum((bsAssociation(:, t) == b)' & opBelongs == v);
        end
    end
end
for b = 1:nBS
    loadDist(:, b) = loadDist(:, b) / sum(loadDist(:, b));
    mask(b) = (sum(loadDist(:, b)) > 0);
end
loadDist(isnan(loadDist)) = 0;
end

