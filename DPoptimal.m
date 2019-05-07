function [userRates, userFraction, btd] = DPoptimal(netSettings, opSettings, ...
    capacityPerUser, bs, minReq, sliceCats)
% DPoptimal Social optimal bidding/rate allocation under DIFFPRICE resource
% sharing criterion.
% Params:
%   NetSettings: network profile
%   OpSettings: operator profile
%   capacityPerUser: capacity perceived per user
%   bs: base station association vector, bs(u) is the BS serving user u.
%   minReq: 1 x V, minimal rate required by slices
%   sliceCats: 1 x V, indicating whether a slice is an elastic slice or
%   not.
% Return:
%   userRates: perceived rate of each user.
%   userFraction: the fraction of time (of associated bs) allocated to each user.
%   btd: perceived user BTDs.

nUsers = netSettings.users;
nSlices = size(opSettings.s_o, 2);    
nBasestations = netSettings.bsNS;
shareVec = opSettings.s_o;
opBelongs = opSettings.ops_belongs;

for v = 1:nSlices
    if(sliceCats(v) == 1)
        minReq(opBelongs == v) = 0;
    end   
end

shareDist = opSettings.shareDist;
userFraction = zeros(1, nUsers);

% Optimize bid based on bidtoutility
options = optimoptions('fmincon','Display','off','Algorithm','interior-point');
% bidding constraint
constMat = zeros(nSlices, nUsers);
constVec = shareVec';
initialBid = 1e-5 * ones(nUsers, 1);

for v = 1:nSlices
    constMat(v, :) = (opBelongs == v)';
    initialBid(opBelongs == v) = shareVec(v) / sum(opBelongs == v);
end

optimBid = fmincon(@(x) -dpbidtoutil(x', bs, opBelongs, capacityPerUser, ...
    shareVec, shareDist, minReq, sliceCats), initialBid, constMat, ...
    constVec, [], [], 1e-5 * ones(nUsers, 1), nBasestations * ones(nUsers, 1), ...
    @(x) minrateconstraint(netSettings, opSettings, x, bs, ...
    capacityPerUser, minReq), options);

cBid = optimBid';

% Compute the rate allocation
for b = 1:nBasestations
    lb = sum(cBid(bs == b));
    if (lb <= 1)
        userFraction(bs == b) = cBid(bs == b) ./ lb;
    else
        for v = 1:nSlices
            lvb = sum(cBid(bs == b & opBelongs == v));
            if (lvb <= shareDist(v, b))
                userFraction(bs == b & opBelongs == v) = cBid(bs == b & opBelongs == v);
            else
                totalSurplus = 1;
                totalWeight = 0;
                for cSlice = 1:nSlices
                    totalWeight = totalWeight + max(0, sum(cBid(opBelongs == cSlice ...
                        & bs == b)) - shareDist(cSlice, b));
                    totalSurplus = totalSurplus - min(shareDist(cSlice, b), ...
                        sum(cBid(opBelongs == cSlice & bs == b)));
                end
                fvb = shareDist(v, b) + totalSurplus / totalWeight * (lvb ...
                    - shareDist(v, b));
                userFraction(bs == b & opBelongs == v) = cBid(bs == b & ...
                    opBelongs == v) .* fvb ./ lvb;
            end
        end
    end
end

userRates = userFraction .* capacityPerUser;
btd=1 ./ userRates;

%---------- Begin nested function --------------------
    function [C, Ceq] = minrateconstraint(netSettings, opSettings, ...
            cBid, bs_, ...
            capacityPerUser, minReq)
    nUsers_ = netSettings.users;
    nSlices_ = size(opSettings.s_o, 2);
    nBasestations_ = netSettings.bsNS;
    shareVec_ = opSettings.s_o;
    opBelongs_ = opSettings.ops_belongs;
    shareDist_ = opSettings.shareDist;
    userFraction_ = zeros(1, nUsers_);
    Ceq = [];
    for b_ = 1:nBasestations_
        lb_ = sum(cBid(bs_ == b_));
        if (lb_ <= 1)
            userFraction_(bs_ == b_) = cBid(bs_ == b_) ./ lb_;
        else
            for v_ = 1:nSlices_
                lvb_ = sum(cBid(bs_ == b_ & opBelongs_ == v_));
                if (lvb_ <= shareDist_(v_, b_))
                    userFraction_(bs_ == b_ & opBelongs_ == v_) = cBid(bs_ == b_ & opBelongs_ == v_);
                else
                    totalSurplus_ = 1;
                    totalWeight_ = 0;
                    for cSlice_ = 1:nSlices_
                        totalWeight_ = totalWeight_ + max(0, sum(cBid(opBelongs_ == cSlice_ ...
                            & bs_ == b_)) - shareDist_(cSlice_, b_));
                        totalSurplus_ = totalSurplus_ - min(shareDist_(cSlice_, b_), ...
                            sum(cBid(opBelongs_ == cSlice_ & bs_ == b_)));
                    end
                    fvb_ = shareDist_(v_, b_) + totalSurplus_ / totalWeight_ * (lvb_ ...
                        - shareDist_(v_, b_));
                    userFraction_(bs_ == b_ & opBelongs_ == v_) = cBid(bs_ == b_ & ...
                        opBelongs_ == v_) .* fvb_ ./ lvb_;
                end
            end
        end
    end
    userRates_ = userFraction_ .* capacityPerUser;
    C = (minReq - userRates_)';
    end
%---------- End nested function --------------------
end