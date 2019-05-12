function [userRates, userFraction, btd, nRounds, isViolation] = DIFFPRICE(netSettings, opSettings, ...
    capacityPerUser, bs, minReq, waterfilling, sliceCats, phi)
% Share constrained sharing with guarantee
% Resources are provisioned according to DIFFPRICE.
% 
% Params:
%   NetSettings: network profile
%   OpSettings: operator profile
%   capacityPerUser: capacity perceived per user
%   bs: base station association vector, bs(u) is the BS serving user u.
%   minReq: 1 x V, minimal rate required by slices
%   waterfilling: logical, whether to use waterfilling or equal allocating
%   surplus.
%   sliceCats: 1 x V, slice types, 1 for non-inelastic, 0 for inelastic.
%   phi: 1 x nUsers, weight for each user.
% Return:
%   userRates: perceived rate of each user.
%   userFraction: the fraction of time (of associated bs) allocated to each user.
%   btd: perceived user BTDs.
%   nRounds: number of rounds needed to converge.
%   isViolation: 1 if the module is > 1, 0 otherwise.

if (nargin == 7)
    phi = ones(1, netSettings.users);
end

nUsers = netSettings.users;
nSlices = size(opSettings.s_o, 2);
nBasestations = netSettings.bsNS;
shareVec = opSettings.s_o;
opBelongs = opSettings.ops_belongs;
shareDist = opSettings.shareDist;
userFraction = zeros(1, nUsers);
eps = 1e-30; % threshold for convergence. 

% Initialize by equal weight. Log utility is not defined at zero.
cBid = zeros(1, nUsers);

for v = 1:nSlices
    cBid(opBelongs == v) = shareVec(v) / sum(opBelongs == v);
end
prevBid = cBid + 1;
cnt = 0;

while(norm(prevBid - cBid) > eps) 
    prevBid = cBid;
    for v = 1:nSlices
        nextBid = diffpriceiteration(cBid, v, shareDist, shareVec, ...
            opBelongs, bs, capacityPerUser, minReq(v), waterfilling, ...
            sliceCats(v), phi);
    end
    if(~all(nextBid > 0))
        assert(all(nextBid > 0), 'Unexpected negative bids.')
    end
    cBid = nextBid;
    cnt = cnt + 1;
    if (cnt == 40)
        disp('break to prevent looping forever, convergence might not achieved.');
    end
end

nRounds = cnt;
% check whether is a violation
rmax = 0;
phimax = 0;
for v = 1:nSlices
    for b = 1:nBasestations
        rmax = max(rmax, sum(minReq(v) ./ capacityPerUser(opBelongs == v ...
            & bs == b)) * sum(opBelongs...
            == v & bs == b));
        phimax = max(phimax, sum(opBelongs == v & bs == b) / sum(opBelongs == v));
    end
end

% compute effective B
effB = 0;
for v = 1:nSlices
    effB = max(effB, length(unique(bs(opBelongs == v))));
end
module = nSlices * (rmax / (1 - rmax)) * (1 + effB * phimax);
isViolation = module > 1 || module < -1;

% Provision resources
for b = 1:nBasestations
    lb = sum(cBid(bs == b));
    if (lb <= 1)
        userFraction(bs == b) = cBid(bs == b) / lb;
    else
        for v = 1:nSlices
            lvb = sum(cBid(bs == b & opBelongs == v));
            if (lvb <= shareDist(v, b))
                userFraction(bs == b & opBelongs == v) = cBid(bs == b ...
                    & opBelongs == v);
            else
                totalSurplus = 1;
                totalWeight = 0;
                for cSlice = 1:nSlices
                    totalWeight = totalWeight + max(0, sum(cBid(opBelongs == cSlice ...
                        & bs == b)) - shareDist(cSlice, b));
                    totalSurplus = totalSurplus - min(shareDist(cSlice, b), ...
                        sum(cBid(opBelongs == cSlice & bs == b)));
                end
                fvb = (shareDist(v, b) + totalSurplus * (lvb - shareDist(v, b)) / totalWeight);
                userFraction(bs == b & opBelongs == v) = fvb ...
                    .* (cBid(bs == b & opBelongs == v)) ...
                    ./ lvb;
                    %./ sum(cBid(bs == b & opBelongs == v));
            end
        end
    end
end

userRates = userFraction .* capacityPerUser;
btd=1 ./ userRates;

end