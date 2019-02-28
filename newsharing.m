function [userRates, userFraction, btd] = newsharing(NetSettings, OpSettings, ...
    capacityPerUser, bs, thresholdFactor)
% New sharing scheme under test. Heuristic: more loaded a BS is, more of
% its resource should be committed to fulfilling the minimal requirements
% of users.
% But first, let's try pure PS.
V = size(OpSettings.s_o, 2);
nBasestations = NetSettings.bsNS;
shareVec = OpSettings.s_o;
weights = nBasestations * OpSettings.w_i;
opBelongs = OpSettings.ops_belongs;
nUsers = NetSettings.users;
baseline = nUsers / nBasestations;
userFraction = zeros(1, nUsers);

% (TODO: a bit of different, each slice should has its own GPS factor, instead 
% of a common one.)
for b = 1:nBasestations
    userIdx = find(bs == b); %index of associated users
    GPSFactor = zeros(1, V);
    % judge based on significant presence
    nvb = zeros(1, V);
    for v = 1:V
        nvb(v) = sum(bs == b & opBelongs == v);
        GPSFactor(v) = nvb(v) > (thresholdFactor * shareVec(v) * baseline);
    end
    ops = unique(OpSettings.ops_belongs(bs == b)); % active tenants
    for u = userIdx
        if (GPSFactor(opBelongs(u)))
            qd = (opBelongs == opBelongs(u)) .* (bs == b);
            userFraction(u) = shareVec(opBelongs(u)) / sum(shareVec(ops)) ...
                / sum(qd); 
        else
            remainingShare = 1;
            for v = ops
                remainingShare = remainingShare - GPSFactor(v) ...
                    * shareVec(v) / sum(shareVec(ops));
            end
            userFraction(u) = remainingShare * OpSettings.w_i(u) ...
                ./ sum(OpSettings.w_i(bs == b));
        end
    end
end
% 
% for u = 1:NetSettings.users
%     nvb = sum(bs == bs(u) & opBelongs == opBelongs(u));
%     GPSFactor = 0;
%     if (nvb > OpSettings.s_o(OpSettings.ops_belongs(u)) * baseline)
%         GPSFactor = 1; % if significant presence, claim its share
%     end
%     ops = unique(OpSettings.ops_belongs([bs==bs(u)])); % active tenants
%     currentO = OpSettings.ops_belongs(u); % tenant of user u
%     qd = ([OpSettings.ops_belongs == currentO] .* [bs==bs(u)]);
%     userFraction(u) = GPSFactor * OpSettings.s_o(OpSettings.ops_belongs(u))...
%         / sum(OpSettings.s_o(ops)) / sum(qd) + (1 - GPSFactor) * ...
%         OpSettings.w_i(u) ./ sum(OpSettings.w_i(bs == bs(u)));
% end

userRates = userFraction .* capacityPerUser;
btd=1 ./ userRates;

end

