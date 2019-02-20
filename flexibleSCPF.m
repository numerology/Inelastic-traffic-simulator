function [userRates, userFraction, btd] = flexibleSCPF(NetSettings, OpSettings, ...
    capacityPerUser, bs)
%flexibleSCPF New version of SCPF that allows slices to specify different share
%allocation across different BSs.
% 
% Params:
%   NetSettings: network profile
%   OpSettings: operator profile
%   capacityPerUser: capacity perceived per user
%   bs: base station association vector, bs(u) is the BS serving user u.
% Return:
%   userRates: perceived rate of each user.
%   userFraction: the fraction of time (of associated bs) allocated to each user.
%   btd: perceived user BTDs.

V = size(OpSettings.s_o, 2);
nBasestations = NetSettings.bsNS;
shareVec = OpSettings.s_o;
shareDist = OpSettings.shareDist; % share distribution, V x B
opBelongs = OpSettings.ops_belongs;

% Compute sum of rescaled load
rescaledLoad = zeros(1, V);
for v = 1:V
    for b = 1:nBasestations
        rescaledLoad(v) = rescaledLoad(v) + shareDist(v, b) * sum(bs(:) == b ...
            & opBelongs(:) == v);
    end
end

% (TODO: optimize this)
for u = 1:NetSettings.users
    curV = opBelongs(u);
    curB = bs(u);
    sumRescaledLoad = 0;
    for v = 1:V
        sumRescaledLoad = sumRescaledLoad + shareVec(v) * shareDist(v, curB) ...
            * sum(bs == curB & opBelongs == v) / rescaledLoad(v);
    end
    userFraction(u) = shareDist(curV, curB) * shareVec(curV) / rescaledLoad(curV)...
        / sumRescaledLoad;
end
userRates = userFraction .* capacityPerUser;
btd=1 ./ userRates;
end

