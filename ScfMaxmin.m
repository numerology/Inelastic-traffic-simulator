% Author: JxZheng
function [rates, fractions, btd] = ScfMaxmin(NetSettings, OpSettings, ...
    capacityPerUser, bs, userDemands)
% Compute the rate allocation under Share Constrained Fair (SCF) with alpha
% = infty, equivalent to weighted maxmin case.
% In static slicing, first all the resources are splitted according to the
% share allocation among slices, then a maxmin is performed to determine
% the time fraction served to each user. Finally we compute the service
% rate according to the time fraction.

% NetSettings: network profile
% OpSettings: operator profile
% c_u: capacity perceived per user (deprecate)
% bs: base station association vector, bs(u) is the BS serving user u.
% userDemands: R x nUsers, userDemands(:, u) is the demand vector of user
% u.
% return:
% rates: vector, rates(u) = rate served to user u
% fraction: fraction of frame/time allocated to user u (deprecate)
% btd: btd perceived by user 

originalCap = [ones(19, 1); 3 * ones(7, 1)]; % total capacities.