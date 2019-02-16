function [loadDist] = getloaddistribution(opSettings, netSettings, bsAssociation, ...
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

V = opSettings.operators;
nBS = netSettings.bsNS;
opBelongs = opSettings.ops_belongs;
loadDist = zeros(V, nBS);
for t = 1:simulationTime
    for v = 1:V
        for b = 1:nBS
            loadDist(v, b) = loadDist(v, b) ...
                + sum((bsAssociation(:, t) == b)' & opBelongs == v);
        end
    end
end
loadDist = loadDist / simulationTime;
for v = 1:V
    loadDist(v, :) = loadDist(v, :) / sum(loadDist(v, :));
end
end

