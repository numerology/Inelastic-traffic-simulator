function [userRates, userFraction, btd] = MAXWEIGHT_BR(netSettings, ...
    opSettings, capacityPerUser, bs)
%MAXWEIGHT_BR SCG by first bidding several rounds for user weights. Bidding
%is carried out according to the best response.
% Params:
%   netSettings: network profile
%   opSettings: operator profile
%   capacityPerUser: capacity perceived per user
%   bs: base station association vector, bs(u) is the BS serving user u.
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
eps = 1e-4; % threshold for convergence. 

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
        cBid = bestresponseiteration(cBid, v, shareDist, shareVec, ...
            opBelongs, bs, capacityPerUser);
    end
    cnt = cnt + 1;
    if (cnt == 20)
        break;
    end
    %assert(cnt < 1000, 'bidding does not converge.'); % Cant wait forever.
end
% fprintf('Required %f rounds to converge.\n', cnt);

% Checking if there is base station overcommitted
% for b = 1:nBasestations
%     if(sum(cBid(bs == b)) > 1.0001)
%         b;
%         sum(cBid(bs == b));
%     end
%     assert(sum(cBid(bs == b)) <= 1.0001, 'Base station is overbooked.');
% end

for u = 1:nUsers
    userFraction(u) = cBid(u) / sum(cBid(bs == bs(u)));
end

userRates = userFraction .* capacityPerUser;
btd=1 ./ userRates;

end