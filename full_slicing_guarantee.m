% Script for SCG simulation
% with more realistic mobility model in addition to poisson.
clc, close all, clear all
%parpool('local', 40);
warning('off','all');
%% Set up
nSlices = 4; % num of slices

sat = 5; % U/B (use only integers...)
simulationTime = 20; % seconds

phiLevels = 1;alphas = [1, 1, 1, 1]; % legacy parameters
warmup = 0;
bsN = 19;
sectors = 5;
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

minRateReq = 1 / (sat) * ones(1, nSlices);
minRateReq(3:4) = 6 * minRateReq(3:4);

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
    outageDP = outageDP + sum(r < (perUserMinRateReq - 1e-3));
    nRoundsVec(t) = nRounds;
    violationVec(t) = isViolation;
    
    [r, f, b] = DPoptimal(NetSettings, OpSettings, capacityPerUser(:,t)', ...
        bs(:,t)', perUserMinRateReq, sliceCats);
    rates_DPoptimal(:, t)=r;
    fractions_DPoptimal(:, t)=f;
    btd_DPoptimal(:, t)=b;
    outageDPoptimal = outageDPoptimal + sum(r < (perUserMinRateReq - 1e-3));
    
    [r, f, b] = SCPF(NetSettings, OpSettings, capacityPerUser(:,t)', bs(:,t)');
    rates_SCPF(:, t)=r;
    fractions_SCPF(:, t)=f;
    btd_SCPF(:, t)=b;
    outageSCPF = outageSCPF + sum(r < (perUserMinRateReq - 1e-3));
    
    tmpOpSettings = OpSettings;
    tmpOpSettings.shareDist = gpsShareDist;
    [r, f, b] = flexibleGPS(NetSettings, tmpOpSettings, capacityPerUser(:,t)', ...
        bs(:,t)', perUserMinRateReq);
    rates_GPS(:, t) = r;
    fractions_GPS(:, t)=f;
    btd_GPS(:, t)=b;
    outageGPS = outageGPS + sum(r < (perUserMinRateReq - 1e-3));
    
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
disp('----------------------------------------------------');
fprintf('DP: P(outage) = %d\n', outageDP / totalNumUsers);
fprintf('DP optimal: P(outage) = %d\n', outageDPoptimal / totalNumUsers);
fprintf('SCPF: P(outage) = %d\n', outageSCPF / totalNumUsers);
fprintf('GPS: P(outage) = %d\n', outageGPS / totalNumUsers);

disp('----------------------------------------------------');
fprintf('DP: BTD = %d\n', nanmean(nanmean(btd_DP)));
fprintf('DPoptimal: BTD = %d\n', nanmean(nanmean(btd_DPoptimal)));
fprintf('SCPF: BTD = %d\n', nanmean(nanmean(btd_SCPF)));
fprintf('GPS: BTD = %d\n', nanmean(nanmean(btd_GPS)));

disp('----------------------------------------------------');
for t = 1:simulationTime % stats
    % For each time instant, only account for the set of users receiving at
    % least min rate req under all benchmarks.
    nUsers = NetSettings.users;
    opVec = OpSettings.ops_belongs;

    goodUsers = (rates_GPS(:, t)' > perUserMinRateReq & rates_SCPF(:, t)' ...
        > perUserMinRateReq & rates_DP(:, t)' > perUserMinRateReq ...
        & rates_DPoptimal(:, t)' > perUserMinRateReq);

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

    tmpRatesDPoptimal = nan(size(rates_DPoptimal(:, t)'));
    tmpRatesDPoptimal(goodUsers) = rates_DPoptimal(goodUsers, t);
    utilDPoptimal(t) = ratetoutil(tmpRatesDPoptimal, shareVec, ...
        opVec, sliceCats, perUserMinRateReq);
end
fprintf('DP: util = %d\n', nanmean(utilDP));
fprintf('DPoptimal: util = %d\n', nanmean(utilDPoptimal));
fprintf('SCPF: util = %d\n', nanmean(utilSCPF));
fprintf('GPS: util = %d\n', nanmean(utilGPS));

%% Plot some metrics

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
