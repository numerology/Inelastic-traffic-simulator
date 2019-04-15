function [util] = dpbidtoutil(cBid, bs, opBelongs, capacityPerUser, ...
    shareVec, shareDist)
% dpbidtoutil Given current bid per user, return the social welfare.
%
% Params:
% cBid: current bid per user, 1 x nUsers
% v: current tenant
% bs: base station association vector, 1 x nUsers
% opBelongs: operator association vector, 1 x nUsers
% capacityPerUser: user perceived capacities, 1 x nUsers

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

assert(all(ratesPerUser > 0), 'negative user rates');

util = 0;
for v = 1:nSlices
    util = util + shareVec(v) / sum(opBelongs == v) * sum(log(ratesPerUser(opBelongs == v)));
end

end