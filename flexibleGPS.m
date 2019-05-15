function [userRates, userFraction, btd] = flexibleGPS(netSettings, ...
    opSettings, capacityPerUser, bs, minRateReq)
% GPS resource allocation enabling specifying share allocation per base
% station.
% Params:
%   netSettings, opSettings: simulation setup. opSettings.shareDist is the
%   share allocation across BSs per slice.
%   capacityPerUser: 1 x U, capacity perceived per user
%   bs: 1 x U, base station association vector, bs(u) is the BS serving
%   user u.
%   minRateReq: 1 x U, minimal rate requirement across slices.

nUsers = netSettings.users;
nSlices = size(opSettings.s_o, 2);
nBasestations = netSettings.bsNS;
%shareVec = opSettings.s_o;
opBelongs = opSettings.ops_belongs;
shareDist = opSettings.shareDist;
userFraction = zeros(1, nUsers);

for b = 1:nBasestations
    ops = unique(opBelongs(bs == b));
    sumWeight = sum(shareDist(ops, b));
    for v = 1:nSlices
        nLocalUsers = sum(opBelongs == v & bs == b);
        minRequiredFraction = sum(minRateReq(opBelongs == v & bs == b) ...
            ./ capacityPerUser(opBelongs == v & bs == b));
        minRequiredBid = minRequiredFraction * sumWeight;
        if (minRequiredBid > shareDist(v, b))
            disp('GPS cannot fulfill.')
            % Allocate bid propto minRateReq
            userFraction(opBelongs == v & bs == b) = shareDist(v, b) / sumWeight .* ...
                minRateReq(opBelongs == v & bs == b) ./ capacityPerUser(opBelongs == v & bs == b) ...
                ./ minRequiredFraction;
        else
            % First allocate the minimal, then equal allocation
            surPlus = shareDist(v, b) - minRequiredBid;
            userFraction(opBelongs == v & bs == b) = minRateReq(opBelongs == v & bs == b) ...
                ./ capacityPerUser(opBelongs == v & bs == b) ...
                + surPlus / nLocalUsers / sumWeight; 
        end
    end
end

userRates = userFraction .* capacityPerUser;
btd = 1 ./ userRates;
end