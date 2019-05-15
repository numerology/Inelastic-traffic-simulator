function [util] = ratetoutil(rates, shareVec, opBelongs, sliceCats, ...
    minRateReq, phi, eps)
% rateutil Compute utility function based on rate allocation, and share
% distribution across users.
% Also need to take into account the nan rate allocation.
% Per user utility is truncated from above at eps < 0.
if(nargin == 5)
    phi = ones(size(opBelongs)); % uniform phi by default.
    eps = -100; % uniform eps by default.
end

util = 0;
nUsers = length(opBelongs);
for u = 1:nUsers
    if (sliceCats(opBelongs(u)) == 1)
        util = util + max(shareVec(opBelongs(u)) * phi(u) / sum(phi(opBelongs == opBelongs(u))) ...
            * log(rates(u) - minRateReq(u)), eps);
    end
end
end