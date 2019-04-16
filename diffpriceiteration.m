function [nextBid] = diffpriceiteration(cBid, v, shareDist, shareVec, ...
    opBelongs, bs, capacityPerUser, minRate, waterfilling)
% diffpriceiteration Generate the bid per user when it is slice v's turn to
% bid.
% Params:
%   cBid: nUsers x 1, current bid
%   v: idx of slice performing this iteration
%   shareDist: V x B, share distribution
%   opBelongs: nUsers x 1, slice association per user
%   bs: nUsers x 1, bs association per user
%   capacityPerUser: 1 x nUsers
%   minRate: minimal rate requirement, scalar
%   waterfilling: logical, 1 when want to user waterfilling to allocate
%   bids, otherwise surplus bids will be allocated equally.

nSlices = size(shareDist, 1);
nBasestations = size(shareDist, 2);
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
        minBidReq(b) = aob(b) / (capacity / minRate / nob - 1);
    elseif (nob * minRate / capacity <= shareDist(v, b))
        minBidReq(b) = nob * minRate / capacity;
    else
        totalSurplus = 1 - shareDist(v, b);
        totalWeight = 0;
        for cSlice = 1:nSlices
            if (cSlice == v)
                continue
            end
            totalWeight = totalWeight + max(0, sum(cBid(opBelongs == cSlice ...
                & bs == b)) - shareDist(cSlice, b));
            totalSurplus = totalSurplus - min(shareDist(cSlice, b), ...
                sum(cBid(opBelongs == cSlice & bs == b)));
        end
        minBidReq(b) = shareDist(v, b) + totalWeight / (totalSurplus / ...
            (nob * minRate / capacity - shareDist(v, b)) - 1);
    end
end

% distribute lob to each user.
% assert(shareVec(v) >= sum(minBidReq), 'Insufficient share to support minRate');
if (waterfilling)
    sliceUserDist = zeros(1, nBasestations);
    for b = 1:nBasestations
        sliceUserDist(b) = sum(opBelongs == v & bs == b);
    end
    surplus = shareVec(v) - sum(minBidReq);
    if (surplus >= 0)
        % Can fulfill all the min rate req.
        addition = waterfill(minBidReq, surplus, sliceUserDist); % TODO: weights may vary?
        loVec = minBidReq + addition;
    else
        % Cannot fulfill the minimal rate req.
        % Assign bid propto minBidReq
        loVec = shareVec(v) .* minBidReq ./ sum(minBidReq);
    end
    nextBid = cBid;
    for b = 1:nBasestations
        nextBid(opBelongs == v & bs == b) = loVec(b) / sliceUserDist(b);
    end
    assert(all(nextBid > 0), 'Unexpected negative bid.');
else
    surplusShare = shareVec(v) - sum(minBidReq);
    nextBid = cBid;
    if (surplusShare < 0 || any(minBidReq < 0))
        nextBid(opBelongs == v) = shareVec(v) / sum(opBelongs == v);
    else
        for b = 1:nBasestations
            nextBid(opBelongs == v & bs == b) = minBidReq(b) / sum(opBelongs == v ...
                & bs == b);
        end
        nextBid(opBelongs == v) = nextBid(opBelongs == v) + surplusShare ...
            / sum(opBelongs == v);
    end   
end

end
