function [poutageGain, utilityGain, poutageErrors, utilityErrors, ...
    elasticShareVec] = GREETsimulation(meanFactorVec, ...
    simulationTime, profile, mode)
%GREETsimulation simulate the gain over GPS in utility and gain over SCPF
%in P(outage) under various load, provided the network configuration
%pattern
%   Params:
%   satVector: vector of the multiplicator of mean rate over min rate.
%   simulationTime: duration of the simulation, < 1000
%   profile: index of configuration, 1 to 5. They are:
%   (for 1 - 5, Slices 1 and 2 are rate-adaptive, 3 and 4 are best-effort.)
%   1: all slices are uniform (RWP). Namely: Uniform.
%   2: 1 and 2 SLAW, 3 and 4 RWP. Namely: TBD.
%   3: all slices are SLAW, with the same hotspots. Namely: Heterogeneous
%      aligned.
%   4: all slices are SLAW, with different hotspots. Namely: Heterogeneous
%      orthogonal.
%   5: all slices are SLAW, with different hotspots, also, Slice 1 and 2
%      have different phi's within each slice. Namely: TBD
%   6: change Slice 1 (or 1 and 2) in setting 2 to purely inelastic slice
%      mode: 
%        change Slice 1 for simulation for the P(outage) - utility 
%        trade-off figure; 
%        change Slice 1 and 2 for simulation for the multi-scenario figure.
%      Namely: Mixed
%
%   Returns:
%   
%   elaticShareVec: 1 x length(meanFactorVec), corresponding share of
%   elastic slices.
%   Then dependes on the mode, when mode = 1:
%   poutageGain: 4 x length(satVector), poutage achieved by 4 benchmarks
%   under different load
%   utilityGain: 4 x length(satVector), utility achieved by 4 benchmarks
%   under different load.
%   errors: 4 x length(meanFactorVec), 1.96 std confidence interval (95%
%   confidence).
%   When mode = 2:
%   poutageGain: 1 x length(satVector), poutage under SCPF - poutage under 
%   DP-practical.
%   utilityGain: 1 x length(satVector), utility under DP-practical -
%   utility under GPS.
elasticShareVec = zeros(1, length(meanFactorVec));
poutageErrors = zeros(1, length(meanFactorVec));
utilityErrors = zeros(1, length(meanFactorVec));
if (mode == 1)
    poutageGain = zeros(4, length(meanFactorVec)); 
    utilityGain = zeros(4, length(meanFactorVec)); 
else
    poutageGain = zeros(1, length(meanFactorVec)); % over SCPF
    utilityGain = zeros(1, length(meanFactorVec)); % over GPS
end

nSlices = 4; % num of slices
phiLevels = 1;alphas = [1, 1, 1, 1]; % legacy parameters
warmup = 0;
bsN = 19;
sectors = 3;
interdistance = 20;
outageTol = 0.01;
minSharePerBS = 0.001;
% User mobility patterns:
% RWP for roughly uniform spatial loads.
model = {'RWP'};

for i = 1:length(meanFactorVec)
    sat = 5;
    
    fprintf('Begin ratio %d \n', meanFactorVec(i));
    
    shareVec = 1/4 * ones(1, 4); % this only means user numbers are the same.
    sliceCats = [0 0 1 1];
    
    % Generate network profile.
    if profile <= 5
        adjustedProfile = profile;
    else
        adjustedProfile = 2;
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
    

    minRateReq(3:4) = meanFactorVec(i) * minRateReq(3:4);
    
    
    [shareDist, gpsShareDist, shareVec] = sharedimension(minRateReq, loadDist, outageTol, ...
            minSharePerBS, 1, 0, sliceCats, bsMask, meanCapacityDist);
    elasticShareVec(i) = sum(sliceCats .* shareVec);
    fprintf('Elastic share = %d', elasticShareVec(i));
    
    weightPerUser = zeros(1, NetSettings.users);
    if (profile ~= 5)
        for v = 1:nSlices
            weightPerUser(OpSettings.ops_belongs == v) = shareVec(v) ...
                ./ sum(OpSettings.ops_belongs == v);
        end
        OpSettings.w_i = weightPerUser;
    else
        phi = ones(1, NetSettings.users);
        for v = 1:nSlices
            userSet = find(OpSettings.ops_belongs == v);
            nUserCurSlice = sum(OpSettings.ops_belongs == v);
            cnt = 1;
            % half of the users are given twice weight
            for u = userSet
                phi(u) = 2;
                if(cnt > nUserCurSlice / 2)
                    break
                end
                cnt = cnt + 1;
            end
            weightPerUser(OpSettings.ops_belongs == v) = shareVec(v) ...
                .* phi(OpSettings.ops_belongs == v) ./ ...
                sum(phi(OpSettings.ops_belongs == v));
        end
        OpSettings.w_i = weightPerUser;
    end

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
    outageDP = zeros(1, simulationTime);
    outageSCPF = zeros(1, simulationTime);
    outageGPS = zeros(1, simulationTime);
    outageDPoptimal = zeros(1, simulationTime);
    % initialization
    rates_DP = zeros(NetSettings.users, simulationTime);
    rates_DPoptimal = zeros(NetSettings.users, simulationTime);
    rates_SCPF = zeros(NetSettings.users, simulationTime);
    rates_GPS = zeros(NetSettings.users, simulationTime);
    parfor t=1:simulationTime
        totalNumUsers = totalNumUsers + NetSettings.users;
        if (mode == 1)
            [r, f, b] = DPoptimal(NetSettings, OpSettings, capacityPerUser(:,t)', ...
                bs(:,t)', perUserMinRateReq, sliceCats);
            rates_DPoptimal(:, t)=r;
            outageDPoptimal(t) = sum(r < (perUserMinRateReq - 1e-3)) / NetSettings.users;
        end
        [r, f, b, nRounds, isViolation] = DIFFPRICE(NetSettings, OpSettings, capacityPerUser(:,t)', ...
            bs(:,t)', minRateReq, 0, sliceCats, weightPerUser);
        rates_DP(:, t)=r;
        outageDP(t) = sum(r < (perUserMinRateReq - 1e-5)) / NetSettings.users;

        [r, f, b] = SCPF(NetSettings, OpSettings, capacityPerUser(:,t)', bs(:,t)');
        rates_SCPF(:, t)=r;
        outageSCPF(t) = sum(r < (perUserMinRateReq)) / NetSettings.users;

        tmpOpSettings = OpSettings;
        tmpOpSettings.shareDist = gpsShareDist;
        [r, f, b] = flexibleGPS(NetSettings, tmpOpSettings, capacityPerUser(:,t)', ...
            bs(:,t)', perUserMinRateReq); % dummy minreq.
        rates_GPS(:, t) = r;
        outageGPS(t) = sum(r < (perUserMinRateReq)) / NetSettings.users;

        fprintf('finish at time %d\n', t);
    end
    if (mode == 1)
        poutageGain(1, i) = mean(outageDP);
        poutageGain(2, i) = mean(outageDPoptimal);
        poutageGain(3, i) = mean(outageSCPF);
        poutageGain(4, i) = mean(outageGPS);
    else
        if (outageSCPF == 0)
            poutageGain(i) = 1;
            poutageErrors(i) = 0;
        else
            fprintf('SCPF outage is %d', mean(outageSCPF));
            fprintf('GREET outage is %d', mean(outageDP));
            poutageGain(i) = mean(outageSCPF) / mean(outageDP);
            poutageErrors(i) = 1.96 * (std(outageSCPF) + std(outageDP)) ...
                / sqrt(simulationTime);
        end
    end
    
    
    % compute utility
    for t = 1:simulationTime % stats
        % For each time instant, only account for the set of users receiving at
        % least min rate req under all benchmarks.
        opVec = OpSettings.ops_belongs;
        
        if (mode == 1)
            goodUsers = (rates_GPS(:, t)' > perUserMinRateReq & rates_SCPF(:, t)' ...
                > perUserMinRateReq & rates_DP(:, t)' > perUserMinRateReq ...
                & rates_DPoptimal(:, t)' > perUserMinRateReq);
            
            tmpRatesGPS = nan(size(rates_GPS(:, t)'));
            tmpRatesGPS(goodUsers) = rates_GPS(goodUsers, t);
            utilGPS(t) = ratetoutil_old(tmpRatesGPS, shareVec, ...
                opVec, sliceCats, perUserMinRateReq, weightPerUser);

            tmpRatesSCPF = nan(size(rates_SCPF(:, t)'));
            tmpRatesSCPF(goodUsers) = rates_SCPF(goodUsers, t);
            utilSCPF(t) = ratetoutil_old(tmpRatesSCPF, shareVec, ...
                opVec, sliceCats, perUserMinRateReq, weightPerUser);

            tmpRatesDP = nan(size(rates_DP(:, t)'));
            tmpRatesDP(goodUsers) = rates_DP(goodUsers, t);
            utilDP(t) = ratetoutil_old(tmpRatesDP, shareVec, ...
                opVec, sliceCats, perUserMinRateReq, weightPerUser);
            
            tmpRatesDPoptimal = nan(size(rates_DPoptimal(:, t)'));
            tmpRatesDPoptimal(goodUsers) = rates_DPoptimal(goodUsers, t);
            utilDPoptimal(t) = ratetoutil_old(tmpRatesDPoptimal, shareVec, ...
                opVec, sliceCats, perUserMinRateReq, weightPerUser);
        else
            goodUsers = (rates_GPS(:, t)' > perUserMinRateReq ...
                & rates_DP(:, t)' > perUserMinRateReq);
            tmpRatesGPS = nan(size(rates_GPS(:, t)'));
            tmpRatesGPS(goodUsers) = rates_GPS(goodUsers, t);
            utilGPS(t) = ratetoutil_old(tmpRatesGPS, shareVec, ...
                opVec, sliceCats, perUserMinRateReq, weightPerUser);

            tmpRatesDP = nan(size(rates_DP(:, t)'));
            tmpRatesDP(goodUsers) = rates_DP(goodUsers, t);
            utilDP(t) = ratetoutil_old(tmpRatesDP, shareVec, ...
                opVec, sliceCats, perUserMinRateReq, weightPerUser);
        end
    end
    if (mode == 1)
        utilityGain(1, i) = nanmean(utilDP);
        utilityGain(2, i) = nanmean(utilDPoptimal);
        utilityGain(3, i) = nanmean(utilSCPF);
        utilityGain(4, i) = nanmean(utilGPS);
    else
        utilityGain(i) = nanmean(utilDP) - nanmean(utilGPS);
        utilityErrors(i) = 1.96 * (nanstd(utilDP) + nanstd(utilGPS)) / sqrt(simulationTime);
    end
end

end