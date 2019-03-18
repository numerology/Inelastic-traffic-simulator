function [userRates, userFraction, btd] = BESTRESPONSE(netSettings, ...
    opSettings, capacityPerUser, bs)
% Sharing criterion where each tenant seeks to maximize their own utility
% function unilaterally, which is defined by a weighted sum of log rate.
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
shareDist = OpSettings.shareDist;
eps = 1e-5; % threshold for convergence. 


end

