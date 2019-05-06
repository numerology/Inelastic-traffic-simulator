function [shareDist, gpsShareDist, shareVec] = sharedimension(minRateReq, ...
     meanLoadDist, eps, minShare, prob, binomial, sliceCats, bsMask)
% sharedimension Find a share allocation for each slice such that the
% minimal rate requirement is satisfied with probability > 1-eps.
% Also, for the extra part of the 
% Params:
% minRateReq: 1 x V
% meanRateReq: 1 x V
% meanLoadDist: V x B
% shareVec: 1 x V
% eps: outage tolerance
% minShare: minimal share allocated to a slice at each resource.
% binomial: logical, 1 when dimension according to binomial distribution, 0
% when poisson dist.
% sliceCats: 1 x V, slice type, 0 for inelastic, 1 for elastic.
% bsMask: 1 x B, 0 means the corresponding BS will be excluded from
% computation.

% Return:
% shareDist: V x B
% gpsShareDist: V x B
% shareVec: 1 x V

if (nargin == 7)
    bsMask = ones(1, size(meanLoadDist, 2));
end

V = length(minRateReq);
shareVec = zeros(1, V);
B = size(meanLoadDist, 2);
minimalShare = zeros(size(meanLoadDist));
gpsShareDist = zeros(size(meanLoadDist));
shareDist = zeros(size(meanLoadDist));
for v = 1:V
    for b = 1:B
        if (~bsMask)
            % if not an active BS, every one get 0.
            minimalShare(v, b) = 0;
            continue;
        end
        if (sliceCats(v) == 0)
            if (binomial)
                minimalShare(v, b) = binoinv(1 - eps, meanLoadDist(v, b), prob) ...
                    * minRateReq(v);
            else
                minimalShare(v, b) = poissinv(1 - eps, meanLoadDist(v, b)) ...
                    * minRateReq(v);
            end
            minimalShare(v, b) = max(minShare, minimalShare(v, b));
        else
            minimalShare(v, b) = minShare;
        end
    end
end

for b = 1:B
    assert(sum(minimalShare(:, b)) <= 1, ...
            'Base station is overbooked.');
    if (~bsMask)
        continue;
    end
    for v = 1:V
        if (sliceCats(v) == 0)
            gpsShareDist(v, b) = minimalShare(v, b);
            shareDist(v, b) = minimalShare(v, b);
        else
            if (sum(meanLoadDist(sliceCats == 1, b)) == 0)
                gpsShareDist(v, b) = (1 - sum(minimalShare(sliceCats ...
                    == 0, b))) / sum(sliceCats);
            else
                gpsShareDist(v, b) = (1 - sum(minimalShare(sliceCats == 0, b))) ...
                    * meanLoadDist(v, b) / sum(meanLoadDist(sliceCats == 1, b));
            end
            shareDist(v, b) = minimalShare(v, b);
        end
    end
end

% assert(all(all(gpsShareDist(:, bsMask > 0) > 0)), 'Non positive share allocated.');
for v = 1:V
    shareVec(v) = sum(gpsShareDist(v, bsMask > 0));
end
end