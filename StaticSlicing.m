function [userRates, userFraction, btd] = StaticSlicing(NetSettings, OpSettings, ...
    capacityPerUser, bs, userDemands)
% Compute the rate allocation under static slicing.
% In static slicing, first all the resources are splitted according to the
% share allocation among slices, then a maxmin is performed to determine
% the time fraction served to each user. Finally we compute the service
% rate according to the time fraction.

% NetSettings: network profile
% OpSettings: operator profile
% capacityPerUser: capacity perceived per user
% bs: base station association vector, bs(u) is the BS serving user u.
% userDemands: R x nUsers, userDemands(:, u) is the demand vector of user
% u.
% return:
% rates: vector, rates(u) = rate served to user u
% fraction: fraction of frame/time allocated to user u
% btd: btd perceived by user

originalCap = [ones(19, 1); 3 * ones(7, 1)]; % total capacities.
userFraction = zeros(size(OpSettings.opsBelongs));
for o = 1:OpSettings.operators
    userWeights = (OpSettings.opsBelongs == o);
    % split capacities of front base stations and backhauls
    [cFraction] = maxmin(userWeights', originalCap ...
        * OpSettings.shareVec(o), userDemands);
    userFraction = userFraction + cFraction;
end
userRates = capacityPerUser .* userFraction; % in future, capacity perUser will be determined by link estimation.

btd=1./userRates;
end
