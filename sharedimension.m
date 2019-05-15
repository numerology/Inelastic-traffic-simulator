function [shareDist, gpsShareDist, shareVec] = sharedimension(minRateReq, ...
     meanLoadDist, eps, minShare, prob, binomial, sliceCats, bsMask, ...
     capacityDist)
% sharedimension Find a share allocation for each slice such that the
% minimal rate requirement is satisfied with probability > 1-eps.
% Also, for the extra part of the 
% Params:
% minRateReq: 1 x V: for inelastic slices, this is interpreted as the
% minimal rate requirement per user; for elastic/best-effort slices, this
% represents the average rate expected per user.
% meanLoadDist: V x B
% shareVec: 1 x V
% eps: outage tolerance
% minShare: minimal share allocated to a slice at each resource.
% binomial: logical, 1 when dimension according to binomial distribution, 0
% when poisson dist.
% sliceCats: 1 x V, slice type, 0 for inelastic, 1 for elastic.
% bsMask: 1 x B, 0 means the corresponding BS will be excluded from
% computation.
% capacityDist: V x B, mean capacity seen by slice v's user at BS b, by default
% uniform unit capacity is assumed

% Return:
% shareDist: V x B
% gpsShareDist: V x B
% shareVec: 1 x V

V = length(minRateReq);
B = size(meanLoadDist, 2);

if (nargin == 7)
    bsMask = ones(1, size(meanLoadDist, 2));
    capacityDist = ones(V, B);
end

shareVec = zeros(1, V);
minimalShare = zeros(size(meanLoadDist));
meanShare = zeros(size(meanLoadDist));
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
                    * minRateReq(v) / capacityDist(v, b);
            else
                minimalShare(v, b) = poissinv(1 - eps, meanLoadDist(v, b)) ...
                    * minRateReq(v) / capacityDist(v, b);
            end
            minimalShare(v, b) = max(minShare, minimalShare(v, b));
        else
            minimalShare(v, b) = 0;
        end
    end
end

for v = 1:V
    if (sliceCats(v) == 1)
        for b = 1:B
            if (~bsMask)
                meanShare(v, b) = 0;
                continue;
            end
            % The min rate req of elastic slices is interpreted as mean
            % rate req.
            meanShare(v, b) = meanLoadDist(v, b) * minRateReq(v) / ...
                capacityDist(v, b);
        end
        
    end
end

% For elastic users under GPS:
% First figure out the e^o for each elastic slice, then water filling
% Between elastic slices, share are partitioned proportional to fraction of
% requirement.

elasticNeed = zeros(1, V);
for v = 1:V
    if (sliceCats(v) == 1)
        elasticNeed(v) = sum(meanShare(v, :));
    end
end
totalElasticNeed = sum(elasticNeed);

assert(sum(sum(minimalShare)) + sum(meanShare(:, b)) <= sum(bsMask), ...
    'System ill-dimensioned, reduce the elastic requirements.');

% iteratively fill the minimal gap
remainderPerBS = 1 - sum(minimalShare, 1);
while(totalElasticNeed > 1e-6)
    [minRemainder, minIdx] = min(remainderPerBS);
    for b = 1:B
        maxDeduct = totalElasticNeed / sum(bsMask & (~isnan(remainderPerBS)));
        if (bsMask)
            remainderPerBS(b) = remainderPerBS(b) - min(minRemainder, maxDeduct);
            totalElasticNeed = totalElasticNeed - min(minRemainder, maxDeduct);
            if (remainderPerBS(b) < 1e-6)
                remainderPerBS(b) = nan;
            end
        end
    end
end

elasticDistGPS = zeros(1, B);
for b = 1:B
    if(isnan(remainderPerBS(b)) || remainderPerBS(b) == 0)
        elasticDistGPS(b) = 1 - sum(minimalShare(:, b));
    else
        elasticDistGPS(b) = 1 - sum(minimalShare(:, b)) - remainderPerBS(b);       
    end
end

% GPS share dimensioning
for b = 1:B
    assert(sum(minimalShare(:, b)) + elasticDistGPS(b) <= 1, ...
        sprintf('base station %b is overbooked.', b));
    %assert(sum(minimalShare(:, b)) <= 1, ...
    %       'Base station is overbooked.');
    if (~bsMask)
        continue;
    end
    for v = 1:V
        if (sliceCats(v) == 0) % inelastic / rate-adaptive
            gpsShareDist(v, b) = minimalShare(v, b);
            shareDist(v, b) = minimalShare(v, b);
        else % elastic / best-effort
            % Proportionally allocate elasticDist
            if (sum(meanLoadDist(sliceCats == 1, b)) == 0)
                gpsShareDist(v, b) = minShare;
            else
                gpsShareDist(v, b) = elasticDistGPS(b) * meanLoadDist(v, b) ...
                    / sum(meanLoadDist(sliceCats == 1, b));
            end
            shareDist(v, b) = minShare;
        end
    end
end

% assert(all(all(gpsShareDist(:, bsMask > 0) > 0)), 'Non positive share allocated.');
for v = 1:V
    shareVec(v) = sum(gpsShareDist(v, bsMask > 0));
end
end