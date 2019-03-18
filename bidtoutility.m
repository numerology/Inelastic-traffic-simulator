function [util] = bidtoutility(cBid, v, bs, opBelongs, capacityPerUser, ...
    shareVec)
% bidtoutility Given current bid per user, return the utility function of tenant
% v.
% Params:
% cBid: current bid per user, 1 x nUsers
% v: current tenant
% bs: base station association vector, 1 x nUsers
% opBelongs: operator association vector, 1 x nUsers
% capacityPerUser: user perceived capacities, 1 x nUsers
nUsers = length(cBid);
util = 0;
for u = 1:nUsers
    if (opBelongs(u) == v)
        util = util + log(capacityPerUser(u) * cBid(u) / sum(cBid(bs == bs(u))));
    end
end
util = util * shareVec(v) / sum(opBelongs == v);
end