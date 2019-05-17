function [util] = ratetoutil_old(rates, shareVec, opBelongs, sliceCats, ...
    minRateReq, phi)
% rateutil Compute utility function based on rate allocation, and share
% distribution across users.
% Also need to take into account the nan rate allocation.
% Per user utility is truncated from above at eps < 0.
if(nargin == 5)
    phi = ones(size(opBelongs)); % uniform phi by default.
end

util = 0;
nSlices = length(shareVec);
for slice = 1:nSlices
    if (sum(opBelongs == slice) == 0)
        continue
    end
    if (sliceCats(slice))
        util = util + shareVec(slice) / sum(phi(~isnan(rates))) * nansum(phi(opBelongs == slice) ./ ...
            sum(phi(opBelongs == slice)) .* ...
            log(rates(opBelongs == slice) - minRateReq(opBelongs == slice)));
    end
end
end

