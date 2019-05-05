% Script for SCG simulation
% with more realistic mobility model in addition to poisson.
clc, close all, clear all
parpool('local', 40);
warning('off','all');
%% Set up
nSlices = 4; % num of slices
sat = 3; % U/B (use only integers...)
simulationTime = 1000; % seconds
phiLevels = 1;alphas = [1, 1, 1, 1]; % legacy parameters
warmup = 0;
bsN = 19;
sectors = 3;
interdistance = 1000;
outageTol = 0.05;
minSharePerBS = 0.01;
% User mobility patterns:
% RWP for roughly uniform spatial loads.
model = {'RWP'};
shareVec = 1/4 * ones(1, 4); % this only means user numbers are the same.

%% Mobility and Link estimation
[NetSettings, OpSettings, capacityPerUser, bs, userPos, bsPos] = ...
    networkconfiguration(simulationTime, ...
    warmup, bsN, sectors,...
    interdistance, model,...
    shareVec, phiLevels, sat, o, alphas);

%% Adjust share distribution for new proposed scheme according to the load distribution
% the sum of share across BSs <= share * |B| per slice.
loadDist = getloaddistribution(OpSettings, NetSettings, bs, simulationTime);
% use a similar heuristic to allocate shares

minRateReq = 0.025 * capacity / (varFactor * perBSLoad) * ones(1, nSlice);
minRateReq(3:4) = 0;
[shareDist, gpsShareDist, shareVec] = sharedimension(minRateReq, loadDist, outageTol, ...
        minSharePerBS, varFactor, 0, sliceCats);
    
weightPerUser = zeros(1, NetSettings.users);
for v = 1:nSlice
    weightPerUser(tmpOpSettings.ops_belongs == v) = shareVec(v) ...
        ./ sum(tmpOpSettings.ops_belongs == v);
end
OpSettings.w_i = weightPerUser;

OpSettings.shareDist = loadDist;
OpSettings.s_o = shareVec;

%% Compute fractions
ppm = ParforProgMon('Simulating resource sharing : ', NetSettings.simulation_time);
parfor t=1:simulationTime
    [r, f, b] = DIFFPRICE(NetSettings, OpSettings, capacityPerUser(:,t)', ...
        bs(:,t)', minRateReq, 0);
    [r, f, b] = DPoptimal(NetSettings, OpSettings, capacityPerUser(:,t)', ...
        bs(:,t)', perUserMinRateReq, sliceCats);
    [r, f, b] = SCPF(NetSettings, OpSettings, capacityPerUser(:,t)', bs(:,t)');
    rates_SCPF(:, t)=r;
    fractions_SCPF(:, t)=f;
    btd_SCPF(:, t)=b;
    tmpOpSettings = OpSettings
    tmpOpSettings.shareDist = gpsShareDist;
    [r, f, b] = flexibleGPS(NetSettings, tmpOpSettings, capacityPerUser(:,t)', ...
        bs(:,t)', ones(1, nUsers)); % dummy minreq.
    ratesGPS(:, t) = r;
    ppm.increment();
end