function [C, Ceq] = bidstrata_nonlinear(weight, shareDist, opBelongs, bs, ...
    nSlices, nBasestations, nUsers)
% bidstrata_nonlinear Function used in bidstrata as the nonlinear
% constraints.
dimMat = zeros(nSlices * nBasestations, nUsers);
compDimMat = zeros(nSlices * nBasestations, nUsers);
for v = 1:nSlices
    for b = 1:nBasestations
        dimMat(b * (v - 1) + b, opBelongs == v & bs == b) = 1;
        compDimMat(b * (v - 1) + b, opBelongs ~= v & bs == b) = 1;
    end
end

Ceq = 0;
sobVec = zeros(nSlices * nBasestations, 1);
for v = 1:nSlices
    for b = 1:nBasestations
        sobVec(b * (v - 1) + b) = shareDist(v, b);
    end
end

C = dimMat * weight - max(sobVec, 1 - compDimMat * weight);

end