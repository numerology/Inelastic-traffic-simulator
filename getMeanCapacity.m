function [meanCapacityDist] = getMeanCapacity(opSettings, netSettings,...
    bs, capacityPerUser, simulationTime)
%getMeanCapacity Get the mean capacity of slice v's user at Bs b
%   Params:
% opSettings, netSettings: simulation profiles;
% bs: bs(u, t) is the BS associated with u at time t.
% simulationTime: duration of simulation
V = opSettings.operators;
B = netSettings.bsNS;
meanCapacityDist = zeros(V, B);
opBelongs = repmat(opSettings.ops_belongs', [1 simulationTime]);

for v = 1:V
    for b = 1:B
        meanCapacityDist(v, b) = 1 / nanmean(1 ./ capacityPerUser(opBelongs == v ...
            & bs == b));
        if(isnan(meanCapacityDist(v, b)))
            meanCapacityDist(v, b) = 1; % dummy.
        end
    end
end

end

