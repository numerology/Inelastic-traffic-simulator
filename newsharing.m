function [userRates, userFraction, btd] = newsharing(NetSettings, OpSettings, ...
    capacityPerUser, bs)
% New sharing scheme under test. Heuristic: more loaded a BS is, more of
% its resource should be committed to fulfilling the minimal requirements
% of users.
% But first, let's try pure PS.
V = size(OpSettings.s_o, 2);
nBasestations = NetSettings.bsNS;
shareDist = OpSettings.shareDist; % share distribution, V x B
weights = nBasestations * OpSettings.w_i;
opBelongs = OpSettings.ops_belongs;
nUsers = NetSettings.users;
baseline = nUsers / nBasestations;

for u = 1:NetSettings.users
    nb = sum(bs == bs(u));
    GPSFactor = exp(nb / baseline - 1) / (exp(nb / baseline - 1) + 1);
    ops = unique(OpSettings.ops_belongs([bs==bs(u)])); % active tenants
    currentO = OpSettings.ops_belongs(u); % tenant of user u
    qd = ([OpSettings.ops_belongs == currentO] .* [bs==bs(u)]);
    userFraction(u) = GPSFactor * OpSettings.s_o(OpSettings.ops_belongs(u))...
        / sum(OpSettings.s_o(ops)) / sum(qd) + (1 - GPSFactor) * ...
        OpSettings.w_i(u) ./ sum(OpSettings.w_i(bs == bs(u)));
end

userRates = userFraction .* capacityPerUser;
btd=1 ./ userRates;

end

