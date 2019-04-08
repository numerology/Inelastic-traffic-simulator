function [shareDist] = sharedimension(minRateReq, meanLoadDist, shareVec, ...
    eps, minShare)
% sharedimension Find a share allocation for each slice such that the
% minimal rate requirement is satisfied with probability > 1-eps.
% Params:
% minRateReq: 1 x V
% meanLoadDist: V x B
% shareVec: 1 x V
% minShare: minimal share allocated to a slice at each resource.
% Return:
% shareDist: V x B

V = length(shareVec);
B = size(meanLoadDist, 2);
minimalShare = zeros(size(meanLoadDist));
spareShare = zeros(B, 1);
shareDist = zeros(size(meanLoadDist));
for v = 1:V
    for b = 1:B
        minimalShare(v, b) = poissinv(1 - eps, meanLoadDist(v, b)) ...
            * minRateReq(v);
    end
end
% check 
for v = 1:V
    assert(sum(minimalShare(v, :)) <= shareVec(v), ...
        'Not enough share for slices.');
end
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

end