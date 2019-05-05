% Script for SCG simulation
% with more realistic mobility model in addition to poisson.
clc, close all, clear all
%parpool('local', 40);
warning('off','all');
%% Set up
nSlices = 4; % num of slices
sat = 3; % U/B (use only integers...)
simulationTime = 10; % seconds
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
sliceCats = [0 0 1 1];

%% Mobility and Link estimation
[NetSettings, OpSettings, capacityPerUser, bs, userPos, bsPos] = ...
    networkconfiguration(simulationTime, ...
    warmup, bsN, sectors,...
    interdistance, model,...
    shareVec, phiLevels, sat, nSlices, alphas);

%% Adjust share distribution for new proposed scheme according to the load distribution
% the sum of share across BSs <= share * |B| per slice.
loadDist = getloaddistribution(OpSettings, NetSettings, bs, simulationTime);
% use a similar heuristic to allocate shares

minRateReq = 0.001 * 1 / (sat) * ones(1, nSlices);
minRateReq(3:4) = 0;
[shareDist, gpsShareDist, shareVec] = sharedimension(minRateReq, loadDist, outageTol, ...
        minSharePerBS, 1, 0, sliceCats);
    
weightPerUser = zeros(1, NetSettings.users);
for v = 1:nSlice
    weightPerUser(OpSettings.ops_belongs == v) = shareVec(v) ...
        ./ sum(OpSettings.ops_belongs == v);
end
OpSettings.w_i = weightPerUser;

OpSettings.shareDist = loadDist;
OpSettings.s_o = shareVec;

%% Compute fractions
%ppm = ParforProgMon('Simulating resource sharing : ', NetSettings.simulation_time);
for t=1:simulationTime
    [r, f, b] = DIFFPRICE(NetSettings, OpSettings, capacityPerUser(:,t)', ...
        bs(:,t)', minRateReq, 0);
    rates_DP(:, t)=r;
    fractions_DP(:, t)=f;
    btd_DP(:, t)=b;
    [r, f, b] = DPoptimal(NetSettings, OpSettings, capacityPerUser(:,t)', ...
        bs(:,t)', perUserMinRateReq, sliceCats);
    rates_DPoptimal(:, t)=r;
    fractions_DPoptimal(:, t)=f;
    btd_DPoptimal(:, t)=b;
    [r, f, b] = SCPF(NetSettings, OpSettings, capacityPerUser(:,t)', bs(:,t)');
    rates_SCPF(:, t)=r;
    fractions_SCPF(:, t)=f;
    btd_SCPF(:, t)=b;
    tmpOpSettings = OpSettings
    tmpOpSettings.shareDist = gpsShareDist;
    [r, f, b] = flexibleGPS(NetSettings, tmpOpSettings, capacityPerUser(:,t)', ...
        bs(:,t)', ones(1, nUsers)); % dummy minreq.
    rates_GPS(:, t) = r;
    fractions_GPS(:, t)=f;
    btd_GPS(:, t)=b;
end