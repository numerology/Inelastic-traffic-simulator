function [nextBid] = diffpriceiteration(cBid, v, shareDist, shareVec, ...
    opBelongs, bs, capacityPerUser, minRate)
% diffpriceiteration Generate the bid per user when it is slice v's turn to
% bid.
% Params:
%   cBid: nUsers x 1, current bid
%   v: idx of slice performing this iteration
%   shareDist: V x B, share distribution
%   opBelongs: nUsers x 1, slice association per user
%   bs: nUsers x 1, bs association per user

nSlices = size(shareDist, 1);
nBasestations = size(shareDist, 2);
nUsers = length(opBelongs);
nUsersOfV = sum(opBelongs == v);
capacity = unique(capacityPerUser);
assert(length(capacity) == 1, 'Only allow uniform capacity');

if (nUsersOfV == 0) % in case there is no user on slice v.
    nextBid = cBid;
    return
end

aob = zeros(1, nBasestations);
for slice = 1:nSlices
    if (slice == v)
        continue
    end
    for b = 1:nBasestations
        aob(b) = aob(b) + sum(cBid(opBelongs == slice & bs == b));
    end
end
lobVec = zeros(1, nBasestations);
minBidReq = zeros(1, nBasestations); % minimal bid at each BSs to meet minRate.
for b = 1:nBasestations
    nob = sum(opBelongs == v & bs == b);
    % incrementally test 3 cases to figure out how much bid is needed.
    if (aob(b) / (capacity / minRate / nob - 1) <= 1 - aob(b))
        lobVec(b) = aob(b) / (capacity / minRate / nob - 1);
    elseif (nob * minRate / capacity <= shareDist(v, b))
        lobVec(b) = nob * minRate / capacity;
    else
        totalSurplus = 1 - shareDist(v, b);
        totalWeight = 0;
        for cSlice = 1:nSlices
            totalWeight = totalWeight + max(0, sum(cBid(opBelongs == cSlice & bs == b)) - shareDist(cSlice, b));
            
        end
    end
end
