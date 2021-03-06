function [util] = dpbidtoutil(cBid, bs, opBelongs, capacityPerUser, ...
    shareVec, shareDist, minRateReq, sliceCats)
% dpbidtoutil Given current bid per user, return the social welfare.
%
% Params:
% cBid: current bid per user, 1 x nUsers
% v: current tenant
% bs: base station association vector, 1 x nUsers
% opBelongs: operator association vector, 1 x nUsers
% capacityPerUser: user perceived capacities, 1 x nUsers
% minRateReq: users' minimal rate requirement, 1 x nUsers
% sliceCats: category of slices' services, 1 for elastic user and 2 for
% inelastic user, 1 x nSlices.

nUsers = length(bs);
nSlices = length(shareVec);
nBasestations = size(shareDist, 2);
ratesPerUser = zeros(1, nUsers);

assert(length(cBid) == nUsers, 'Invalid length of bids.');
assert(length(opBelongs) == nUsers, 'Invalid length of opBelongs.');
assert(length(capacityPerUser) == nUsers, 'Invalid length of capacity.');
assert(all(size(shareDist) == [nSlices, nBasestations]));

for b = 1:nBasestations
    lb = sum(cBid(bs == b));
    if (lb <= 1)
        ratesPerUser(bs == b) = capacityPerUser(bs == b) .* cBid(bs == b) ./ lb;
    else
        for v = 1:nSlices
            lvb = sum(cBid(bs == b & opBelongs == v));
            if (lvb <= shareDist(v, b))
                ratesPerUser(bs == b & opBelongs == v) = capacityPerUser(bs == b & opBelongs == v) ...
                    .* cBid(bs == b & opBelongs == v);
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
                ratesPerUser(bs == b & opBelongs == v) = capacityPerUser(bs == b & opBelongs == v) ...
                    .* cBid(bs == b & opBelongs == v) .* fvb ./ lvb;
            end
        end
    end
end

assert(all(ratesPerUser >= 0), 'negative user rates');

% Here to account for minimal rate requirement, one has two options:
% 1. set that in the constraint set and guarantee that the initial point is
% feasible. For example, begin with equal surplus practical approach.
% 2. Penalize a huge amount of utility for the case with user violating the minimal
% rate requirement, for example, -1e5.
% Here for simplicity we use 1, because 2 is too slow.

util = 0;
for v = 1:nSlices
    if (sum(opBelongs == v) == 0)
        continue
    end
    if (sliceCats(v) == 1)
        util = util + shareVec(v) / sum(opBelongs == v) ...
            * nansum(log(ratesPerUser(opBelongs == v & ratesPerUser > minRateReq) ...
            - minRateReq(opBelongs == v & ratesPerUser > minRateReq)));
        util = util - 10000 * sum(ratesPerUser(opBelongs == v) ...
            < minRateReq(opBelongs == v));
    else
        util = util - 10000 * sum(ratesPerUser(opBelongs == v) ...
            < minRateReq(opBelongs == v));
    end
end
end