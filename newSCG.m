function [userRates, userFraction, btd] = newSCG(NetSettings, OpSettings, ...
    capacityPerUser, bs)
%new SCG New version of sharing proposed by Albert.
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
opBelongs = OpSettings.ops_belongs;
remainingPerBs = ones(1, nBasestations);
bidPerUser = zeros(size(OpSettings.w_i));
for v = 1:V
    bidPerUser(opBelongs == v) = shareVec(v) * nBasestations / sum(opBelongs == v);
    lbRate =  shareVec(v) * nBasestations / sum(opBelongs == v) / 2; 
    % Currently, uniform rate requirement, which equals to half of the initial.
    for b = 1:nBasestations
        assert(remainingPerBs(b) > 0, 'Insufficient budget at bs.');
        nvb = sum(bs == b & opBelongs == v);
        if nvb == 0
            continue
        end
        
        if lbRate * nvb > remainingPerBs(b)
            cvb = (lbRate * nvb - remainingPerBs(b)) / nvb;
        else
            cvb = 0;
        end
        
        bidPerUser(bs == b & opBelongs == v) = lbRate - cvb;
        remainingPerBs(b) = remainingPerBs(b) - nvb * (lbRate - cvb);
    end
end

for u = 1:NetSettings.users
    if(bidPerUser(u) == 0)
        userFraction(u) = 1e-6;
    else
        userFraction(u) = bidPerUser(u) / sum(bidPerUser(bs == bs(u) & opBelongs ...
            == opBelongs(u)));
    end
end

userRates = userFraction .* capacityPerUser;
btd=1 ./ userRates;

end