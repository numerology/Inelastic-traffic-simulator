function [inverseP] = inversewaterfill(total, ceiling, nUsers)
% inversewaterfill Inverse water filling function. Try to solve the minmax
% problem of the amount of water in each bucket, subject to the constraint that
% the height of each bucket is upperbounded by vector ceiling.

L = length(ceiling); % number of buckets
totalVol = sum(ceiling);
assert(total <= totalVol, 'No feasible allocation found for water fill.');

ceilPerUser = ceiling ./ nUsers;
[upsortLevels, perm] = sort(ceilPerUser);
nCumUsers = cumsum(nUsers(perm));
nCumUsers = nCumUsers(end) - [0 nCumUsers(1:(L-1))];
[a, b] = size(upsortLevels);
if (a > b)
    upsortLevels = upsortLevels';
    ceilPerUser = ceilPerUser';
    nCumUsers = nCumUsers';
end
extendedLevels = [0 upsortLevels];
delta = [0 cumsum((extendedLevels(2:L) - extendedLevels(1:(L-1))) .* ...
    nCumUsers(1:(L - 1)))];
l = sum(total >= delta); % how many ceils can be reached + 1
level = (total - delta(l)) / nCumUsers(l) + extendedLevels(l);

inverseP = ceiling;
for l = 1:L
    if(level < ceilPerUser(l))
        inverseP(l) = level * nUsers(l);
    end
end
end