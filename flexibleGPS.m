function [userRates, userFraction, btd] = flexibleGPS(netSettings, ...
    opSettings, capacityPerUser, bs)
% GPS resource allocation enabling specifying share allocation per base
% station.
% Params:
%   netSettings, opSettings: simulation setup. opSettings.shareDist is the
%   share allocation across BSs per slice.
%   capacityPerUser: 1 x U, capacity perceived per user
%   bs: 1 x U, base station association vector, bs(u) is the BS serving
%   user u.

for u=1:netSettings.users
    b = bs(u);
    ops = unique(opSettings.ops_belongs([bs==bs(u)]));
    slice = opSettings.ops_belongs(u);
    % number of users on the same slice at the same bs.
    qd = sum(opSettings.ops_belongs == slice & bs == bs(u)); 
    assert(opSettings.shareDist(slice, b) > 0, 'zero allocation for active user');
    userFraction(u) = opSettings.shareDist(slice, b)/...
                                sum(opSettings.shareDist(ops, b))/...
                                sum(qd);
end
userRates = userFraction .* capacityPerUser;
btd = 1 ./ userRates;
end