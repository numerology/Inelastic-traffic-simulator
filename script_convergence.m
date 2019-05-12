% script for convergence examination
% Script for SCG simulation
% with more realistic mobility model in addition to poisson.
clc, close all, clear all
parpool('local', 40);
warning('off','all');
%% Set up
nSlices = 4; % num of slices

sat = 1; % U/B (use only integers...)
simulationTime = 100; % seconds

phiLevels = 1;alphas = [1, 1, 1, 1]; % legacy parameters
warmup = 0;
bsN = 4;
sectors = 3;
interdistance = 20;
outageTol = 0.01;
minSharePerBS = 0.001;
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
    shareVec, phiLevels, sat, nSlices, alphas, 2);

%% Adjust share distribution for new proposed scheme according to the load distribution
% the sum of share across BSs <= share * |B| per slice.
[loadDist, bsMask] = getloaddistribution(OpSettings, NetSettings, bs, ...
    simulationTime);
meanCapacityDist = getMeanCapacity(OpSettings, NetSettings, bs, capacityPerUser, ...
    simulationTime);
% use a similar heuristic to allocate shares

minRateReq = 5 / (sat) * ones(1, nSlices);
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
    
    fprintf('finish at time %d\n', t);
end

%% Plot cdf of convergence
datestring = datestr(now, 30);

figure(2)
hold on; grid on
cdf1 = cdfplot(nRoundsVec(violationVec == 0));
cdf2 = cdfplot(nRoundsVec(violationVec == 1));
set(cdf1, 'marker', '+');
set(cdf2, 'marker', 'o');
legend('When module < 1', 'When module > 1');
xlabel('number of rounds');
ylabel('CDF');
savefig(sprintf('figs/cdf-convergence-%s.fig', datestring));


%% 
