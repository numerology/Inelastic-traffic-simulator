function [userRates, userFraction, btd] = DPoptimal(netSettings, opSettings, ...
    capacityPerUser, bs, minReq, sliceCats)
% DPoptimal Social optimal bidding/rate allocation under DIFFPRICE resource
% sharing criterion.
% Params:
%   NetSettings: network profile
%   OpSettings: operator profile
%   capacityPerUser: capacity perceived per user
%   bs: base station association vector, bs(u) is the BS serving user u.
%   minReq: 1 x V, minimal rate required by slices
%   sliceCats: 1 x V, indicating whether a slice is an elastic slice or
%   not.
% Return:
%   userRates: perceived rate of each user.
%   userFraction: the fraction of time (of associated bs) allocated to each user.
%   btd: perceived user BTDs.

nUsers = netSettings.users;
nSlices = size(opSettings.s_o, 2);
nBasestations = netSettings.bsNS;
shareVec = opSettings.s_o;
opBelongs = opSettings.ops_belongs;
shareDist = opSettings.shareDist;
userFraction = zeros(1, nUsers);

% Optimize bid based on bidtoutility
options = optimoptions('fmincon','Display','off','Algorithm','sqp');
% bidding constraint
constMat = zeros(nSlices, nUsers);
constVec = shareVec';
for v = 1:nSlices
    constMat(v, :) = (opBelongs == v)';
end
initialBid = 1e-5 * ones(nUsers, 1);
optimBid = fmincon(@(x) -dpbidtoutil(x', bs, opBelongs, capacityPerUser, ...
    shareVec, shareDist, minReq, sliceCats), initialBid, constMat, ...
    constVec, [], [], 1e-5 * ones(nUsers, 1), ones(nUsers, 1), [], options);

cBid = optimBid';

% Compute the rate allocation
for b = 1:nBasestations
    lb = sum(cBid(bs == b));
    if (lb <= 1)
        userFraction(bs == b) = cBid(bs == b) ./ lb;
    else
        for v = 1:nSlices
            lvb = sum(cBid(bs == b & opBelongs == v));
            if (lvb <= shareDist(v, b))
                userFraction(bs == b & opBelongs == v) = cBid(bs == b & opBelongs == v);
            else
                totalSurplus = 1;
                totalWeight = 0;
                for cSlice = 1:nSlices
                    totalWeight = totalWeight + max(0, sum(cBid(opBelongs == cSlice ...
                        & bs == b)) - shareDist(cSlice, b));
                    totalSurplus = totalSurplus - min(shareDist(cSlice, b), ...
                        sum(cBid(opBelongs == cSlice & bs == b)));
                end
                fvb = shareDist(v, b) + totalSurplus / totalWeight * (lvb ...
                    - shareDist(v, b));
                userFraction(bs == b & opBelongs == v) = cBid(bs == b & ...
                    opBelongs == v) .* fvb ./ lvb;
            end
        end
    end
end

userRates = userFraction .* capacityPerUser;
btd=1 ./ userRates;

end