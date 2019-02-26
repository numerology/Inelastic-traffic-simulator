function [userRates, userFraction, btd] = newSCG(NetSettings, OpSettings, ...
    capacityPerUser, bs)
%new SCG New version of sharing proposed by Albert.
% 
% Params:
%   NetSettings: network profile
%   OpSettings: operator profile
%   capacityPerUser: capacity perceived per user
%   bs: base station association vector, bs(u) is the BS serving user u.
% Return:
%   userRates: perceived rate of each user.
%   userFraction: the fraction of time (of associated bs) allocated to each user.
%   btd: perceived user BTDs.

V = size(OpSettings.s_o, 2);
nBasestations = NetSettings.bsNS;
shareVec = OpSettings.s_o;
opBelongs = OpSettings.ops_belongs;
remainingPerBs = ones(1, nBasestations);
bidPerUser = zeros(size(OpSettings.w_i));
shareDist = OpSettings.shareDist;

for v = 1:V
    bidPerUser(opBelongs == v) = shareVec(v) * nBasestations...
        / sum(opBelongs == v);
    % Update remaining budget per BS.
    for b = 1:nBasestations
        remainingPerBs(b) = 1 - sum(bidPerUser(opBelongs ~= v & bs == b));
    end
    % Currently, uniform rate requirement.
    for b = 1:nBasestations
        %assert(remainingPerBs(b) > 0, 'Insufficient budget at bs.');
        nvb = sum(bs == b & opBelongs == v);
        if nvb == 0
            continue
        end
        
        if (sum(bidPerUser(opBelongs == v & bs == b)) > ...
                max(shareDist(v, b), remainingPerBs(b)))
            overFlow = sum(bidPerUser(opBelongs == v & bs == b)) - ...
                max(shareDist(v, b), remainingPerBs(b));
            bidPerUser(opBelongs == v & bs == b) = ...
                bidPerUser(opBelongs == v & bs == b) - overFlow / nvb;
            % When redistribute the bid, we can only touch the later BSs,
            % because we don't want to violate the constraints for former
            % BSs.
            nFollowing = sum(bs > b & opBelongs == v);
            bidPerUser(bs > b & opBelongs == v) = ...
                bidPerUser(bs > b & opBelongs == v) + overFlow / nFollowing;
        end
    end
end

for v = 1:V
    % Update remaining budget per BS.
    for b = 1:nBasestations
        remainingPerBs(b) = 1 - sum(bidPerUser(opBelongs ~= v & bs == b));
    end
    % Currently, uniform rate requirement.
    for b = nBasestations:1
        %assert(remainingPerBs(b) > 0, 'Insufficient budget at bs.');
        nvb = sum(bs == b & opBelongs == v);
        if nvb == 0
            continue
        end
        
        if (sum(bidPerUser(opBelongs == v & bs == b)) > ...
                max(shareDist(v, b), remainingPerBs(b)))
            overFlow = sum(bidPerUser(opBelongs == v & bs == b)) - ...
                max(shareDist(v, b), remainingPerBs(b));
            bidPerUser(opBelongs == v & bs == b) = ...
                bidPerUser(opBelongs == v & bs == b) - overFlow / nvb;
            % When redistribute the bid, we can only touch the later BSs,
            % because we don't want to violate the constraints for former
            % BSs.
            nFollowing = sum(bs < b & opBelongs == v);
            bidPerUser(bs < b & opBelongs == v) = ...
                bidPerUser(bs < b & opBelongs == v) + overFlow / nFollowing;
        end
    end
end

for u = 1:NetSettings.users
    if(bidPerUser(u) == 0)
        userFraction(u) = 1e-6; % Prevent infinity
    else
        userFraction(u) = bidPerUser(u) / sum(bidPerUser(bs == bs(u)));
    end
end

userRates = userFraction .* capacityPerUser;
btd=1 ./ userRates;

end