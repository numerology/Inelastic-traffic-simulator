function [userRates, userFraction, btd] = SCG(NetSettings, OpSettings, ...
    capacityPerUser, bs)
% Share constrained sharing with guarantee
% Currently only works under parallel resource usage. The the userDemands will
% be overrided.
% 
% Params:
% NetSettings: network profile
% OpSettings: operator profile
% capacityPerUser: capacity perceived per user
% bs: base station association vector, bs(u) is the BS serving user u.
% Return:
% userRates: perceived rate of each user.
% userFraction: the fraction of time (of associated bs) allocated to each user.
% btd: perceived user BTDs.

assert(isequal(NetSettings.model, {'RWP'}), 'Only random waypoint model is supported.');
V = size(OpSettings.s_o, 2);
nBasestations = NetSettings.bsNS;
shareDist = OpSettings.shareDist; % share distribution, V x B
weights = nBasestations * OpSettings.w_i;
opBelongs = OpSettings.ops_belongs;

% (TODO:optimize the computation here.)
for u = 1:NetSettings.users
    lb = sum(weights(bs(:)==bs(u)));
    if (lb <= 1)
        userFraction(u) = weights(u) ./ lb;
    else
        %disp('lb greater than 1');
        slice = opBelongs(u);
        lvb = sum(weights(bs(:) == bs(u) & opBelongs(:) == slice));
        if (lvb <= shareDist(slice, bs(u)))
            userFraction(u) = weights(u);
        else
            surplus = 1;
            totalMargin = 0;
            for v = 1:V
                localLvb = sum(weights(bs(:) == bs(u) & opBelongs(:) == v));
                surplus = surplus - min(localLvb, shareDist(v, bs(u)));
                totalMargin = totalMargin + max(0, localLvb ...
                    - shareDist(v, bs(u)));
            end
            userFraction(u) = weights(u) / lvb * (shareDist(slice, bs(u)) ...
                + (lvb - shareDist(slice, bs(u))) / totalMargin * surplus);
        end
    end
end

userRates = userFraction .* capacityPerUser;
btd=1 ./ userRates;

end