function [ rates ] = maxminsolve(mmweights, mmcapacities, demands)
% Solve the maxmin rate allocation problem under the resource constraint.
% 
% max  min_i r_i/w_i
% s.t. sum_i d(i, j) r_i <= c_j for each j
%
% Assuming everyone has nonzero weights and demands.
% --------------------------------
% Parameters:
% --------------------------------
% - weights: array of user weights, nUsers x 1 matrix
% - capacities: array of resources capacities, R x 1 matrix
% - demands: resource demands per user, R x nUsers matrix
% --------------------------------
% Ret:
% --------------------------------
% - rate: rate allocation to each user, 1 x nUsers matrix

N = size(mmweights, 2); % number of users.
if N == 0
    rates = [0];
    return
end

rates = zeros(1, N);
R = size(mmcapacities, 1);
remainingCap = mmcapacities;
bottlenecked = zeros(1,N); % none of the users are bottlenecked
bottleneckedResource = zeros(1,R);

while(max(remainingCap) > 1e-7 && min(bottlenecked) == 0 && ...
        min(bottleneckedResource) == 0)
    activeRemainingCap = remainingCap;
    activeRemainingCap(bottleneckedResource > 0) = inf;
    weightedDemands = zeros(R, 1);
    for i = 1:N
        for r = 1:R
            % if any of the resource it needs is bottlenecked,
            % this user is also BNed
            if (bottleneckedResource(r) > 0 && demands(r, i) > 0)
                bottlenecked(i) = 1;
            end
        end
        if (bottlenecked(i) > 0)
            continue;
        end
        for r = 1:R
            weightedDemands(r) = weightedDemands(r) + mmweights(i) ...
                * demands(r, i);
        end
    end
    [increment, minidx] = min(activeRemainingCap ./ weightedDemands);
    if(increment == inf) % to deal with the case demand orthogonal to 
                         % remaining capacity
        break
    end
    bottleneckedResource(minidx) = 1;
    for i = 1:N
        if(bottlenecked(i) > 0)
            continue;
        end
        rates(i) = rates(i) + increment * mmweights(i);
        % all the users have demand in minidx are bned
        if(demands(minidx, i) > 0)
            bottlenecked(i) = 1;
        end
    end

    for r = 1:R
        if(remainingCap == 0)
            continue;
        end
        remainingCap(r) = remainingCap(r) - weightedDemands(r) * ...
            increment;
        if(abs(remainingCap(r)) < 1e-10)
            remainingCap(r) = 0;
            bottleneckedResource(r) = 1; 
            % needed if two resources get BNed at the same time.
        end
    end      
end

end
