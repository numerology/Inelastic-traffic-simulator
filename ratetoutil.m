function [util] = ratetoutil(rates, shareVec, opBelongs, sliceCats, ...
    minRateReq)
% rateutil Compute utility function based on rate allocation, and share
% distribution across users.
util = 0;
nSlice = length(shareVec);
for slice = 1:nSlice
    if (sliceCats(v))
        nCurUsers = sum(opBelongs == slice);
        util = util + shareVec(slice) / nCurUsers * sum(...
            log(rates(opBelongs == slice) - minRateReq(opBelongs == slice)));
    end
end
end