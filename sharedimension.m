function [shareDist] = sharedimension(minRateReq, meanLoadDist, shareVec, ...
    eps, minShare, prob, binomial, GPS)
% sharedimension Find a share allocation for each slice such that the
% minimal rate requirement is satisfied with probability > 1-eps.
% Params:
% minRateReq: 1 x V
% meanLoadDist: V x B
% shareVec: 1 x V
% minShare: minimal share allocated to a slice at each resource.
% binomial: logical, 1 when dimension according to binomial distribution, 0
% when poisson dist.
% GPS: logical, 1 when performing share dimensioning for GPS, i.e., will
% adopt the convention that s^e = 0.
% Return:
% shareDist: V x B

V = length(shareVec);
B = size(meanLoadDist, 2);
minimalShare = zeros(size(meanLoadDist));
spareShare = zeros(B, 1);
shareDist = zeros(size(meanLoadDist));
for v = 1:V
    for b = 1:B
        if (binomial)
            minimalShare(v, b) = binoinv(1 - eps, meanLoadDist(v, b), prob) ...
                * minRateReq(v);
        else
            minimalShare(v, b) = poissinv(1 - eps, meanLoadDist(v, b)) ...
                * minRateReq(v);
        end

        
    end
end
% check 
for v = 1:V
    assert(sum(minimalShare(v, :)) <= shareVec(v), ...
        'Not enough share for slices.');
end

if (GPS)
    for b = 1:B
        assert(sum(minimalShare(:, b)) + V * minShare <= 1, ...
            'Base station is overbooked.');
        spareShare(b) = 1 - sum(minimalShare(:, b)) - V * minShare;
        % currently, allocate spare share prop to minimal allocation.
        for v = 1:V
            if (sum(minimalShare(:, b)) == 0) % if no minimal share is needed, equal allocation.
                shareDist(v, b) = spareShare(b) / V + minShare;
            else
                shareDist(v, b) = minimalShare(v, b) + minShare + spareShare(b) ...
                    * minimalShare(v, b) / sum(minimalShare(:, b));
            end
        end
    end
else
    % if not GPS, only allocate the min share needed. 
    for b = 1:B
        assert(sum(minimalShare(:, b)) + V * minShare <= 1, ...
            'Base station is overbooked.');
        for v = 1:V
            if (sum(minimalShare(:, b)) == 0) % if no minimal share is needed, equal allocation.
                shareDist(v, b) = spareShare(b) / V + minShare;
            else
                shareDist(v, b) = minimalShare(v, b);
            end
        end
    end
end


end