function [bidPerUser] = bidstrata(netSettings, opSettings, capacityPerUser, ...
    bs, minReq)
% bidstrata Figure out the bid for each user given the current system
% profile. The bidding strategy is driven by a maxmin problem seeking to
% maximize the minimal margin above the rate requirement across users.

% Params:
%   netSettings, opSettings: system profile
%   capacityPerUser: 1 x nUsers, capacity perceived by each users.
%   bs: bs(u) is the idx of the base station user u is associated with.
%   minReq: minReq(u) is the minimal rate requirement of user u.






end