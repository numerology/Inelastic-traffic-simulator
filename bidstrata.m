function [bidPerUser] = bidstrata(netSettings, opSettings, capacityPerUser, ...
    bs, minReq)
% bidstrata Figure out the bid for each user given the current system
% profile. The bidding strategy is driven by a maxmin problem seeking to
% maximize the minimal margin above the rate requirement across users.

% Params:
%   netSettings, opSettings: system profile
%   capacityPerUser: 1 x nUsers, capacity perceived by each users.
%   bs: bs(u) is the idx of the base station user u is associated with.
%   minReq: minReq(u) is the minimal rate requirement of user u.
nUsers = size(capacityPerUser, 2);
shareVec = opSettings.s_o;
nSlices = size(shareVec, 2);
opBelongs = OpSettings.ops_belongs;
nBasestations = netSettings.bsNS;

marginFun = @(weight) minReq ./ capacityPerUser - weight;
% Formulate the constraints
% per slice share constraint
shareMat = zeros(nSlices, nUsers);
for v = 1:nSlices
    shareMat(v, opBelongs == v) = 1;
end

% well-dimension constraint at each slice
% slice v at basestation b will be coded at b * (v - 1) + b
dimMat = zeros(nSlices * nBasestations, nUsers);
for v = 1:nSlices
    for b = 1:nBasestations
        dimMat(b * (v - 1) + b, opBelongs == v & bs == b) = 1;
    end
end

bidPerUser = fminimax(marginFun, zeros(size(capacityPerUser)), shareMat, ...
    shareVec', [], [], zeros(size(capacityPerUser)), ...
    ones(size(capacityPerUser)), @(x) bidstrata_nonlinear(x, shareDist, ...
    opBelongs, bs, nSlices, nBasestations, nUsers));

end