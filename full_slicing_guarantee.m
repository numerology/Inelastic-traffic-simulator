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

minRateReq = 0.2 * 1 / (sat) * ones(1, nSlices);
minRateReq(3:4) = 0;
[shareDist, gpsShareDist, shareVec] = sharedimension(minRateReq, loadDist, outageTol, ...
        minSharePerBS, 1, 0, sliceCats);
    
weightPerUser = zeros(1, NetSettings.users);
for v = 1:nSlices
    weightPerUser(OpSettings.ops_belongs == v) = shareVec(v) ...
        ./ sum(OpSettings.ops_belongs == v);
end
OpSettings.w_i = weightPerUser;

OpSettings.shareDist = loadDist;
OpSettings.s_o = shareVec;

perUserMinRateReq = ones(1, NetSettings.users);
for v  = 1:nSlices
    perUserMinRateReq(OpSettings.ops_belongs == v) = minRateReq(v);
end

%% Compute fractions
%ppm = ParforProgMon('Simulating resource sharing : ', NetSettings.simulation_time);
parfor t=1:simulationTime
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
    tmpOpSettings = OpSettings;
    tmpOpSettings.shareDist = gpsShareDist;
    [r, f, b] = flexibleGPS(NetSettings, tmpOpSettings, capacityPerUser(:,t)', ...
        bs(:,t)', ones(1, NetSettings.users)); % dummy minreq.
    rates_GPS(:, t) = r;
    fractions_GPS(:, t)=f;
    btd_GPS(:, t)=b;
    fprintf('finish at time %d\n', t);
end

%% Plot some metrics
datestring = datestr(now, 30);

benchmarks = {'DIFFPRICE-equal surplus', 'GPS', 'SCPF', 'DIFFPRICE-optimal'};
bmWoGPS = {'DIFFPRICE-equal surplus', 'SCPF', 'DIFFPRICE-optimal'};

i1=2;
i2=30;
i3=134;
figure(1);
subplot(3,1,1)
plot(btd_DP(i1,:),'-.b')
hold on
plot(btd_GPS(i1,:),'--r')
plot(btd_SCPF(i1,:),'-g')
plot(btd_DPoptimal(i1,:),':k')

subplot(3,1,2)
plot(btd_DP(i2,:),'-.b')
hold on
plot(btd_GPS(i2,:),'--r')
plot(btd_SCPF(i2,:),'-g')
plot(btd_DPoptimal(i2,:),':k')
subplot(3,1,3)

plot(btd_DP(i3,:),'-.b')
hold on
plot(btd_GPS(i3,:),'--r')
plot(btd_SCPF(i3,:),'-g')
plot(btd_DPoptimal(i3,:),':k')
legend(benchmarks);
savefig(sprintf('figs/btd-trajectory-%s.fig', datestring));

%% 