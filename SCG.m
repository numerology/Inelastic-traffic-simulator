function [userRates, userFraction, btd] = SCG(NetSettings, OpSettings, ...
    capacityPerUser, bs, userDemands)
% Share constrained sharing with guarantee
% Currently only works under parallel resource usage. The the userDemands will
% be overrided.
%
% NetSettings: network profile
% OpSettings: operator profile
% capacityPerUser: capacity perceived per user
% bs: base station association vector, bs(u) is the BS serving user u.
% userDemands: R x nUsers, userDemands(:, u) is the demand vector of user
% u. (obsolete)

end