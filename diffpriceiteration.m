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

nUsers = length(cBid);
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
minBidReqPerUser = zeros(1, nUsers);
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
% If there is no other slices, possible to have minBidReq == 0.
% Set a lowerbound to prevent this.
minBidReq(minBidReq < 1e-4) = 1e-4;
assert(all(minBidReq > 0), 'Unexpected negative minBidReq.');
for b = 1:nBasestations
    minBidReqPerUser(opBelongs == v & bs == b) = minBidReq(b) ...
        / sum(opBelongs == v & bs == b);
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
        idx = (minBidReqPerUser > 0);
        filteredAddition = waterfill(minBidReqPerUser(idx), surplus, ...
            ones(size(minBidReqPerUser(idx))));
        addition = zeros(1, nUsers);
        addition(idx) = filteredAddition;
        loVec = minBidReqPerUser + addition;
    else
        % Cannot fulfill the minimal rate req.
        % Assign bid propto minBidReq
        loVec = shareVec(v) .* minBidReqPerUser ./ sum(minBidReqPerUser);
    end
    nextBid = cBid;
    nextBid(opBelongs == v) = loVec(opBelongs == v);
    assert(all(nextBid > 0), 'Unexpected negative bid.');
else
    surplusShare = shareVec(v) - sum(minBidReq);
    nextBid = cBid;
    if (surplusShare < 0 || any(minBidReq < 0))
        % Cannot satisfy all minimal requirement
        % Knock out users according to MaxSubSet criterion.
        % Until sob can meet the min rate requirements.
        % nextBid(opBelongs == v) = shareVec(v) / sum(opBelongs == v);
        minBidReqPerUser = nan(1, nUsers);
        for b = 1:nBasestations
            if (aob(b) / (capacity / minRate - 1) <= 1 - aob(b))
                minBidReqPerUser(opBelongs == v & bs == b) = ...
                    aob(b) / (capacity / minRate - 1);
            elseif (minRate / capacity <= shareDist(v, b))
                minBidReqPerUser(opBelongs == v & bs == b) = ...
                    minRate / capacity;
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
                minBidReqPerUser(opBelongs == v & bs == b) ...
                    = shareDist(v, b) + totalWeight / (totalSurplus / ...
                    (minRate / capacity - shareDist(v, b)) - 1);
            end
        end
        cTotalShareNeeded = 0;
        admissionControl = zeros(1, nUsers);
        while(cTotalShareNeeded <= shareVec)
            minBidReqPerUser(admissionControl > 0) = inf;
            [minBid, minIdx] = min(minBidReqPerUser);
            admissionControl(minIdx) = 1;
            nextBid(minIdx) = minBid;
            cTotalShareNeeded = cTotalShareNeeded + minBid;
            % update minbid needed
            cBaseStation = bs(minIdx);
            nob = sum(opBelongs == v & bs == cBaseStation ...
                & admissionControl) + 1; % need to count for itself. 
            
            if (aob(cBaseStation) / (capacity / minRate / nob - 1) ...
                    <= 1 - aob(cBaseStation))
                minBidReqPerUser(opBelongs == v & bs == cBaseStation) ...
                   = aob(cBaseStation) / (capacity / minRate / nob - 1);
            elseif (nob * minRate / capacity <= shareDist(v, cBaseStation))
                minBidReqPerUser(opBelongs == v & bs == cBaseStation) ...
                    = nob * minRate / capacity;
            else
                totalSurplus = 1 - shareDist(v, cBaseStation);
                totalWeight = 0;
                for cSlice = 1:nSlices
                    if (cSlice == v)
                        continue
                    end
                    totalWeight = totalWeight + max(0, sum(cBid(opBelongs == cSlice ...
                        & bs == cBaseStation)) - shareDist(cSlice, cBaseStation));
                    totalSurplus = totalSurplus - min(shareDist(cSlice, cBaseStation), ...
                        sum(cBid(opBelongs == cSlice & bs == cBaseStation)));
                end
                minBidReqPerUser(opBelongs == v & bs == cBaseStation) ...
                    = shareDist(v, cBaseStation) + totalWeight / (totalSurplus / ...
                    (nob * minRate / capacity - shareDist(v, cBaseStation)) - 1);
            end
        end
        
        nextBid(opBelongs == v & ~admissionControl) = 1e-5; % give an epsilon bid.
        
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
