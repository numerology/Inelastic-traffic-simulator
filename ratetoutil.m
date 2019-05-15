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
nSlices = length(shareVec);
nUsers = length(opBelongs);
for u = 1:nUsers
    util = util + max(shareVec * phi(u) / sum(opBelongs == opBelongs(u)) ...
        * log(rates(u) - minRateReq(u)), eps);
end
end