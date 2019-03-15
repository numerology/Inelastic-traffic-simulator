function [userRates, userFraction, btd] = MAXWEIGHT(NetSettings, ...
    OpSettings, capacityPerUser, bs)
%biddingSCG SCG by first bidding for two rounds for user weights
% Params:
%   NetSettings: network profile
%   OpSettings: operator profile
%   capacityPerUser: capacity perceived per user
%   bs: base station association vector, bs(u) is the BS serving user u.
% Return:
%   userRates: perceived rate of each user.
%   userFraction: the fraction of time (of associated bs) allocated to each user.
%   btd: perceived user BTDs.

nUsers = NetSettings.users;
nSlices = size(OpSettings.s_o, 2);
nBasestations = NetSettings.bsNS;
shareVec = OpSettings.s_o;
opBelongs = OpSettings.ops_belongs;
remainingPerBs = ones(1, nBasestations);
bidPerUser = zeros(size(OpSettings.w_i));
shareDist = OpSettings.shareDist;
nRounds = 10; % Number of bidding rounds needed. Empirical value. Subject to nSlices. 

% First do two rounds of bidding. Initialize by either all zero or equal
% weight.
cBid = zeros(1, nUsers); % all zero.

% equal weight.
% for v = 1:nSlices
%     cBid(opBelongs == v) = shareVec(v) / sum(opBelongs == v);
% end

for round = 1:nRounds
    for v = 1:nSlices
        cBid = biditeration(cBid, v, shareDist, shareVec, opBelongs, bs);
    end
end

% Checking if there is base station overcommitted
for b = 1:nBasestations
    if(sum(cBid(bs == b)) > 1.0001)
        b
        sum(cBid(bs == b))
    end
    %assert(sum(cBid(bs == b)) <= 1.0001, 'Base station is overbooked.');
end

for u = 1:nUsers
    userFraction(u) = cBid(u); % / sum(cBid(bs == bs(u)));
end

userRates = userFraction .* capacityPerUser;
btd=1 ./ userRates;

end

