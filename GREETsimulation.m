function [poutageGain, utilityGain] = GREETsimulation(satVector, ...
    simulationTime, profile)
%GREETsimulation simulate the gain over GPS in utility and gain over SCPF
%in P(outage) under various load, provided the network configuration
%pattern
%   Params:
%   satVector: vector of sat := |U|/|B|
%   simulationTime: duration of the simulation, < 1000
%   profile: index of configuration, 1 to 5. They are:
%   for 1 - 5, Slices 1 and 2 are rate-adaptive, 3 and 4 are best-effort.
%   1: all slices are uniform (RWP)
%   2: 1 and 2 SLAW, 3 and 4 RWP
%   3: all slices are SLAW, with the same hotspots
%   4: all slices are SLAW, with different hotspots
%   5: all slices are SLAW, with different hotspots, also, Slice 1 and 2
%   have different phi's within each slice
%   6: change Slice 1 (or 1 and 2) in setting 2 to purely inelastic slice

poutageGain = zeros(1, length(satVector)); % over SCPF
utilityGain = zeros(1, length(satVector)); % over GPS

nSlices = 4; % num of slices
phiLevels = 1;alphas = [1, 1, 1, 1]; % legacy parameters
warmup = 0;
bsN = 19;
sectors = 3;
interdistance = 15;
outageTol = 0.01;
minSharePerBS = 0.001;
% User mobility patterns:
% RWP for roughly uniform spatial loads.
model = {'RWP'};
shareVec = 1/4 * ones(1, 4); % this only means user numbers are the same.

for i = 1:length(satVector)
    sat = satVector(i); % U/B (use only integers...)
    sliceCats = [0 0 1 1];
    
    % Generate network profile.
    if profile <= 5
        adjustedProfile = profile;
    else
        adjustedProfile = 5;
    end
    [NetSettings, OpSettings, capacityPerUser, bs, userPos, bsPos] = ...
        networkconfiguration(simulationTime, ...
        warmup, bsN, sectors,...
        interdistance, model,...
        shareVec, phiLevels, sat, nSlices, alphas, adjustedProfile);
    
    % Adjust share dimensioning
    [loadDist, bsMask] = getloaddistribution(OpSettings, NetSettings, bs, ...
        simulationTime);
    meanCapacityDist = getMeanCapacity(OpSettings, NetSettings, bs, capacityPerUser, ...
        simulationTime);
        minRateReq = 1 / (sat) * ones(1, nSlices);
    minRateReq(3:4) = 5 * minRateReq(3:4);

    [shareDist, gpsShareDist, shareVec] = sharedimension(minRateReq, loadDist, outageTol, ...
            minSharePerBS, 1, 0, sliceCats, bsMask, meanCapacityDist);

    weightPerUser = zeros(1, NetSettings.users);
    for v = 1:nSlices
        weightPerUser(OpSettings.ops_belongs == v) = shareVec(v) ...
            ./ sum(OpSettings.ops_belongs == v);
    end
    OpSettings.w_i = weightPerUser;

    OpSettings.shareDist = shareDist;
    OpSettings.s_o = shareVec;

    perUserMinRateReq = zeros(1, NetSettings.users);
    for v  = 1:nSlices
        if (sliceCats(v) == 0)
            perUserMinRateReq(OpSettings.ops_belongs == v) = minRateReq(v);
        end
    end
    minRateReq(sliceCats > 0) = 0;
    
    if (profile ~= 6)
        sliceCats = [1 1 1 1];
    else
        sliceCats = [0 1 1 1]; % or [0 0 1 1];
    end
    
    totalNumUsers = 0;
    outageDP = 0;
    outageSCPF = 0;
    outageGPS = 0;
    
    parfor t=1:simulationTime
        totalNumUsers = totalNumUsers + NetSettings.users;
        [r, f, b, nRounds, isViolation] = DIFFPRICE(NetSettings, OpSettings, capacityPerUser(:,t)', ...
            bs(:,t)', minRateReq, 0, sliceCats);
        rates_DP(:, t)=r;
        fractions_DP(:, t)=f;
        btd_DP(:, t)=b;
        outageDP = outageDP + sum(r < (perUserMinRateReq - 1e-5));

        [r, f, b] = SCPF(NetSettings, OpSettings, capacityPerUser(:,t)', bs(:,t)');
        rates_SCPF(:, t)=r;
        fractions_SCPF(:, t)=f;
        btd_SCPF(:, t)=b;
        outageSCPF = outageSCPF + sum(r < (perUserMinRateReq));

        tmpOpSettings = OpSettings;
        tmpOpSettings.shareDist = gpsShareDist;
        [r, f, b] = flexibleGPS(NetSettings, tmpOpSettings, capacityPerUser(:,t)', ...
            bs(:,t)', ones(1, NetSettings.users)); % dummy minreq.
        rates_GPS(:, t) = r;
        fractions_GPS(:, t)=f;
        btd_GPS(:, t)=b;
        outageGPS = outageGPS + sum(r < (perUserMinRateReq));

        fprintf('finish at time %d\n', t);
    end
    poutageGain(i) = (outageSCPF - outageDP) / totalNumUsers;
    
    % compute utility
    for t = 1:simulationTime % stats
        % For each time instant, only account for the set of users receiving at
        % least min rate req under all benchmarks.
        opVec = OpSettings.ops_belongs;

        goodUsers = (rates_GPS(:, t)' > perUserMinRateReq & rates_SCPF(:, t)' ...
            > perUserMinRateReq & rates_DP(:, t)' > perUserMinRateReq);

        tmpRatesGPS = nan(size(rates_GPS(:, t)'));
        tmpRatesGPS(goodUsers) = rates_GPS(goodUsers, t);
        utilGPS(t) = ratetoutil(tmpRatesGPS, shareVec, ...
            opVec, sliceCats, perUserMinRateReq);

        tmpRatesSCPF = nan(size(rates_SCPF(:, t)'));
        tmpRatesSCPF(goodUsers) = rates_SCPF(goodUsers, t);
        utilSCPF(t) = ratetoutil(tmpRatesSCPF, shareVec, ...
            opVec, sliceCats, perUserMinRateReq);

        tmpRatesDP = nan(size(rates_DP(:, t)'));
        tmpRatesDP(goodUsers) = rates_DP(goodUsers, t);
        utilDP(t) = ratetoutil(tmpRatesDP, shareVec, ...
            opVec, sliceCats, perUserMinRateReq);
    end
    
    utilityGain(i) = nanmean(utilDP) - nanmean(utilGPS);


end

end