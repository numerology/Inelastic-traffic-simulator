function [outputArg1,outputArg2] = GREETsimulation(satVector, profile)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
outputArg1 = inputArg1;
outputArg2 = inputArg2;


nSlices = 4; % num of slices
simulationTime = 1000; % seconds
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
sliceCats = [0 0 1 1];

for i = 1:length(satVector)
    sat = satVector(i); % U/B (use only integers...)
    
    % Begin all uniform case
    [NetSettings, OpSettings, capacityPerUser, bs, userPos, bsPos] = ...
        networkconfiguration(simulationTime, ...
        warmup, bsN, sectors,...
        interdistance, model,...
        shareVec, phiLevels, sat, nSlices, alphas, 1);
    
end

%% Mobility and Link estimation
[NetSettings, OpSettings, capacityPerUser, bs, userPos, bsPos] = ...
    networkconfiguration(simulationTime, ...
    warmup, bsN, sectors,...
    interdistance, model,...
    shareVec, phiLevels, sat, nSlices, alphas, 2);

%% Adjust share distribution for new proposed scheme according to the load distribution
% the sum of share across BSs <= share * |B| per slice.
[loadDist, bsMask] = getloaddistribution(OpSettings, NetSettings, bs, ...
    simulationTime);
meanCapacityDist = getMeanCapacity(OpSettings, NetSettings, bs, capacityPerUser, ...
    simulationTime);
% use a similar heuristic to allocate shares

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
sliceCats = [1 1 1 1];

%% Compute fractions
%ppm = ParforProgMon('Simulating resource sharing : ', NetSettings.simulation_time);
totalNumUsers = 0;
outageDP = 0;
outageDPoptimal = 0;
outageSCPF = 0;
outageGPS = 0;
nRoundsVec = zeros(1, simulationTime);
violationVec = zeros(1, simulationTime);
parfor t=1:simulationTime
    totalNumUsers = totalNumUsers + NetSettings.users;
    [r, f, b, nRounds, isViolation] = DIFFPRICE(NetSettings, OpSettings, capacityPerUser(:,t)', ...
        bs(:,t)', minRateReq, 0, sliceCats);
    rates_DP(:, t)=r;
    fractions_DP(:, t)=f;
    btd_DP(:, t)=b;
    outageDP = outageDP + sum(r < (perUserMinRateReq - 1e-5));
    nRoundsVec(t) = nRounds;
    violationVec(t) = isViolation;
    
    [r, f, b] = DPoptimal(NetSettings, OpSettings, capacityPerUser(:,t)', ...
        bs(:,t)', perUserMinRateReq, sliceCats);
    rates_DPoptimal(:, t)=r;
    fractions_DPoptimal(:, t)=f;
    btd_DPoptimal(:, t)=b;
    outageDPoptimal = outageDPoptimal + sum(r < (perUserMinRateReq));
    
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

end