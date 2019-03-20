function [nextBid] = bestresponseiteration(cBid, v, shareDist, shareVec, ...
    opBelongs, bs, capacityPerUser)
% bestresponseiteration Generate the bid per user on slice v given current bid, 
% together with the network profile, where each tenant maxmize its utility
% function, (currently weighted sum of log rate) unilaterally.
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

userDist = zeros(1, nBasestations);
for b = 1:nBasestations
    userDist(b) = sum(opBelongs == v & bs == b);
end

% Make the matrix for constraint, which rows are divided as follows:
% nSlices * nBasestations : number of local bid constraint
% nSlices: aggregated bid of a tenant <= its share.
% The linear noneq constraint will be constMat * cBid <= constVec
constMat = zeros(nSlices * nBasestations + nSlices, nUsers);
constVec = zeros(nSlices * nBasestations + nSlices, 1);
idx = 1;
for slice = 1:nSlices
    for basestation = 1:nBasestations
        constMat(idx, opBelongs == slice & bs == basestation) = 1;
        constVec(idx) = max(shareDist(slice, basestation), sum(cBid(...
            opBelongs ~= slice & bs == basestation)));
        idx = idx + 1;
    end
end
for slice = 1:nSlices
    constMat(idx, opBelongs == slice) = 1;
    constVec(idx) = shareVec(slice);
    idx = idx + 1;
end

% The eq constraint will be bid of other tenant does not change.
eqMat = zeros(nUsers - nUsersOfV, nUsers);
eqVec = zeros(nUsers - nUsersOfV, 1);
idx = 1;
for user = 1:nUsers
    if (opBelongs(user) ~= v)
        eqMat(idx, user) = 1;
        eqVec(idx) = cBid(user);
        idx = idx + 1;
    end
end

options = optimoptions('fmincon','Display','off');
initial = 1e-3 * ones(size(cBid))';
nextBid = fmincon(@(x) -bidtoutility(x, v, bs, opBelongs, capacityPerUser, ...
    shareVec), initial, constMat, constVec, eqMat, eqVec, zeros(size(initial)), [], [], ...
    options);
nextBid = nextBid';
end