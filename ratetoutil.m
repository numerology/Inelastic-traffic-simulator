function [util] = ratetoutil(rates, shareVec, opBelongs, sliceCats, ...
    minRateReq, phi)
% rateutil Compute utility function based on rate allocation, and share
% distribution across users.
% Also need to take into account the nan rate allocation.
if(nargin == 5)
    phi = ones(size(opBelongs)); % uniform phi by default.
end

util = 0;
nSlice = length(shareVec);
for slice = 1:nSlice
    if (sum(opBelongs == slice) == 0)
        continue
    end
    if (sliceCats(slice))
        util = util + shareVec(slice) * nansum(phi(opBelongs == slice) ./ ...
            sum(phi(opBelongs == slice)) .* ...
            log(rates(opBelongs == slice) - minRateReq(opBelongs == slice)));
    end
end
end