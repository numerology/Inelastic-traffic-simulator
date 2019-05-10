function [nextBid] = diffpriceiteration(cBid, v, shareDist, shareVec, ...
    opBelongs, bs, capacityPerUser, minRate, waterfilling, sliceType)
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
%   sliceType: 0 for inelastic, 1 for non-inelastic.

nUsers = length(cBid);
nSlices = size(shareDist, 1);
nBasestations = size(shareDist, 2);
nUsersOfV = sum(opBelongs == v);

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

minBidReq = zeros(1, nBasestations); % minimal bid at each BSs to meet minRate.
minBidReqPerUser = zeros(1, nUsers);
for b = 1:nBasestations
    totalFractionReq = sum(minRate ./ capacityPerUser(opBelongs == v & bs == b));
    % incrementally test 3 cases to figure out how much bid is needed.
    if (aob(b) / (1 / totalFractionReq - 1) <= 1 - aob(b))
        minBidReq(b) = aob(b) / (1 / totalFractionReq - 1);
    elseif (totalFractionReq <= shareDist(v, b))
        minBidReq(b) = totalFractionReq;
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
            (totalFractionReq - shareDist(v, b)) - 1);
    end
end
% If there is no other slices, possible to have minBidReq == 0.
% Set a lowerbound to prevent this.
minBidReq(minBidReq < 1e-8) = 1e-8;
assert(all(minBidReq > 0), 'Unexpected negative minBidReq.');

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
        disp('Admission control applied.');
        % Cannot satisfy all minimal requirement
        % Knock out users according to MaxSubSet criterion.
        % Until sob can meet the min rate requirements.
        % nextBid(opBelongs == v) = shareVec(v) / sum(opBelongs == v);
        minBidReqPerUser = nan(1, nUsers);
        fractionClaimedPerBS = zeros(1, nBasestations);
        for b = 1:nBasestations
            userSet = find(opBelongs == v & bs == b);
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
            for u = userSet
                if (aob(b) / (1 / (fractionClaimedPerBS(b) + minRate ...
                        / capacityPerUser(u)) - 1) <= 1 - aob(b))
                    minBidReqPerUser(u) = aob(b) / (1 / (fractionClaimedPerBS(b) + minRate ...
                        / capacityPerUser(u)) - 1);
                elseif ((fractionClaimedPerBS(b) + minRate ...
                        / capacityPerUser(u)) <= shareDist(v, b))
                    minBidReqPerUser(u) = fractionClaimedPerBS(b) ...
                        + minRate / capacityPerUser(u);
                else
                    minBidReqPerUser(u) = shareDist(v, b) + totalWeight / (totalSurplus / ...
                        (fractionClaimedPerBS(b) + minRate / capacityPerUser(u) - shareDist(v, b)) - 1);
                end
            end
        end
        minBidReqPerUser(minBidReqPerUser < 1e-5) = 1e-5; % prevent 0 bid.
        minBidReqPerUser(opBelongs ~= v) = inf;
        cTotalShareNeeded = 0;
        admissionControl = zeros(1, nUsers);
        while(cTotalShareNeeded + min(minBidReqPerUser) <= shareVec(v))
            [minBid, minIdx] = min(minBidReqPerUser);
            admissionControl(minIdx) = 1;
            nextBid(minIdx) = minBid;
            cTotalShareNeeded = cTotalShareNeeded + minBid;
            % update minbid needed
            cBaseStation = bs(minIdx);
            fractionClaimedPerBS(cBaseStation) = fractionClaimedPerBS(cBaseStation) ...
                + minRate / capacityPerUser(minIdx);
            userSet = find(opBelongs == v & bs == cBaseStation);
            for u = userSet
                if (aob(cBaseStation) / (1 / (fractionClaimedPerBS(cBaseStation) + minRate ...
                        / capacityPerUser(u)) - 1) <= 1 - aob(cBaseStation))
                    minBidReqPerUser(u) = aob(cBaseStation) / (1 ...
                        / (fractionClaimedPerBS(cBaseStation) + minRate ...
                        / capacityPerUser(u)) - 1);
                elseif ((fractionClaimedPerBS(cBaseStation) + minRate ...
                        / capacityPerUser(u)) <= shareDist(v, cBaseStation))
                    minBidReqPerUser(u) = fractionClaimedPerBS(cBaseStation) ...
                        + minRate / capacityPerUser(u);
                else
                    minBidReqPerUser(u) = shareDist(v, cBaseStation) ...
                        + totalWeight / (totalSurplus / ...
                        (fractionClaimedPerBS(cBaseStation) + minRate ...
                        / capacityPerUser(u) - shareDist(v, cBaseStation)) - 1);
                end
            end
            minBidReqPerUser(minBidReqPerUser < 1e-5) = 1e-5;
            minBidReqPerUser(admissionControl > 0) = inf;
        end
        
        nextBid(opBelongs == v & ~admissionControl) = 1e-8; % give an epsilon bid.
        
    else
        for b = 1:nBasestations
           % should be inversely prop to capacity
            nextBid(opBelongs == v & bs == b) = minBidReq(b) .* ...
                (1 ./ capacityPerUser(opBelongs == v & bs == b)) ...
                ./ sum(1 ./ capacityPerUser(opBelongs == v & bs == b));
        end
        if (sliceType == 1)
            nextBid(opBelongs == v) = nextBid(opBelongs == v) + surplusShare ...
                /sum(opBelongs == v);
        end
    end   
end

end
