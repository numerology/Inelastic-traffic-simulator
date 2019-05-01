function [util] = ratetoutil(rates, shareVec, opBelongs, sliceCats, ...
    minRateReq)
% rateutil Compute utility function based on rate allocation, and share
% distribution across users.
% Also need to take into account the nan rate allocation.
util = 0;
nSlice = length(shareVec);
for slice = 1:nSlice
    if (sum(opBelongs == slice) == 0)
        continue
    end
    if (sliceCats(slice))
        nCurUsers = sum(opBelongs == slice);
        util = util + shareVec(slice) / nCurUsers * nansum(...
            log(rates(opBelongs == slice) - minRateReq(opBelongs == slice)));
    end
end
end