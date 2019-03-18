function [outputBid] = bestresponseiteration(cBid, v, shareDist, shareVec, ...
    opBelongs, bs)
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

% make the matrix for constraint

end