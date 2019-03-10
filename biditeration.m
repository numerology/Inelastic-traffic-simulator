function [outputBid] = biditeration(cBid, v, shareDist, shareVec, opBelongs, bs)
% biditeration Generate the bid per user on slice v given current bid, together
% with the network profile.
% Params:
%   cBid: nUsers x 1, current bid
%   v: idx of slice performing this iteration
%   shareDist: V x B, share distribution
%   opBelongs: nUsers x 1, slice association per user
%   bs: nUsers x 1, bs association per user

% (Jiaxiao 03/09/2019): Here in the first version we only study the case where
% minimal requirement of users belonging to the same tenant are the same.

% Water filling is executed on per base station basis.
% Establish ceiling buckets
nSlices = size(shareDist, 1);
nBasestations = size(shareDist, 2);
nUsers = length(opBelongs);

aob = zeros(1, nBasestations);
for slice = 1:nSlices
    if (slice == v)
        continue
    end
    for b = 1:nBasestations
        aob = aob + sum(cBid(opBelongs == slice & bs == b));
    end
end

ceilBucket = max(shareDist(v, :), 1 - aob);
perBsBid = inversewaterfill(shareVec(v), ceilBucket);
outputBid = cBid;
for b = 1:nBasestations
    outputBid(bs == b & opBelongs == v) = perBsBid(b) / sum((bs == b ...
        & opBelongs == v));
end

end