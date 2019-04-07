function [userRates, userFraction, btd] = DIFFPRICE(NetSettings, OpSettings, ...
    capacityPerUser, bs, minReq)
% Share constrained sharing with guarantee
% Currently only works under parallel resource usage. The the userDemands will
% be overrided.
% 
% Params:
%   NetSettings: network profile
%   OpSettings: operator profile
%   capacityPerUser: capacity perceived per user
%   bs: base station association vector, bs(u) is the BS serving user u.
%   minReq: 1 x V, minimal rate required by slices
% Return:
%   userRates: perceived rate of each user.
%   userFraction: the fraction of time (of associated bs) allocated to each user.
%   btd: perceived user BTDs.

nUsers = netSettings.users;
nSlices = size(opSettings.s_o, 2);
nBasestations = netSettings.bsNS;
shareVec = opSettings.s_o;
opBelongs = opSettings.ops_belongs;
shareDist = opSettings.shareDist;
userFraction = zeros(1, nUsers);
eps = 1e-6; % threshold for convergence. 

% Initialize by equal weight. Log utility is not defined at zero.
cBid = zeros(1, nUsers);

for v = 1:nSlices
    cBid(opBelongs == v) = shareVec(v) / sum(opBelongs == v);
end
prevBid = cBid + 1;
cnt = 0;

while(norm(prevBid - cBid) > eps) 
    prevBid = cBid;
    for v = 1:nSlices
        cBid = diffpriceiteration(cBid, v, shareDist, shareVec, ...
            opBelongs, bs, capacityPerUser, minReq(v));
    end
    cnt = cnt + 1;
    if (cnt == 40)
        disp('break to prevent looping forever, convergence might not achieved.');
    end
end

% Provision resources
for b = 1:nBasestations
    lb = sum(cBid(bs == b));
    if (lb <= 1)
        userFraction(bs == b) = cBid(bs == b) / lb;
    else
        for v = 1:nSlices
            lvb = sum(cBid(bs == b & opBelongs == v));
            if (lvb <= shareDist(v, b))
                userFraction(bs == b & opBelongs == v) = cBid(bs == b ...
                    & opBelongs == v);
            else
                totalSurplus = 1;
                totalWeight = 0;
                for cSlice = 1:nSlices
                    totalWeight = totalWeight + max(0, sum(cBid(opBelongs == cSlice ...
                        & bs == b)) - shareDist(cSlice, b));
                    totalSurplus = totalSurplus - min(shareDist(cSlice, b), ...
                        sum(cBid(opBelongs == cSlice & bs == b)));
                end
                userFraction(bs == b & opBelongs == v) = (shareDist(v, b) + ...
                    totalSurplus * (lvb - shareDist(v, b)) / totalWeight) / ...
                    sum(opBelongs == v & bs == b);
            end
        end
    end
end

% (TODO:optimize the computation here.)
for u = 1:NetSettings.users
    lb = sum(weights(bs(:)==bs(u)));
    if (lb <= 1)
        userFraction(u) = weights(u) ./ lb;
    else
        %disp('lb greater than 1');
        slice = opBelongs(u);
        lvb = sum(weights(bs(:) == bs(u) & opBelongs(:) == slice));
        if (lvb <= shareDist(slice, bs(u)))
            userFraction(u) = weights(u);
        else
            surplus = 1;
            totalMargin = 0;
            for v = 1:nSlices
                localLvb = sum(weights(bs(:) == bs(u) & opBelongs(:) == v));
                surplus = surplus - min(localLvb, shareDist(v, bs(u)));
                totalMargin = totalMargin + max(0, localLvb ...
                    - shareDist(v, bs(u)));
            end
            userFraction(u) = weights(u) / lvb * (shareDist(slice, bs(u)) ...
                + (lvb - shareDist(slice, bs(u))) / totalMargin * surplus);
        end
    end
end

userRates = userFraction .* capacityPerUser;
btd=1 ./ userRates;

end