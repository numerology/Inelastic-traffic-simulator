% script for convergence examination
% Script for SCG simulation
% with more realistic mobility model in addition to poisson.
clc, close all, clear all
%parpool('local', 40);
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
relativeRhoVec = 1 * [[2 2 6 6];
                              [10 10 6 6];
                              [10 10 6 6];
                              [2 2 6 6]]';
nBaseStations = size(relativeRhoVec, 2);
capacity = 1;
sliceCats = [0 0 1 1];
afterSliceCats = [1 1 1 1];
netSettings = [];
netSettings.bsNS = nBaseStations;
opSettings = [];

%% Adjust share distribution for new proposed scheme according to the load distribution
% use a similar heuristic to allocate shares

minRateReq = 0.04 / (sat) * ones(1, nSlices);
minRateReq(3:4) = 1 * minRateReq(3:4);

[shareDist, gpsShareDist, shareVec] = sharedimension(minRateReq, relativeRhoVec, outageTol, ...
        minSharePerBS, 1, 0, sliceCats, ones(1, 4), ones(4, 4));

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
    loadDist = poissrnd(relativeRhoVec);
    nUsers = sum(sum(loadDist));
    totalNumUsers = totalNumUsers + nUsers;
    bsVec = zeros(1, nUsers);
    opVec = zeros(1, nUsers);
    capVec = capacity * ones(1, nUsers);
    ptr = 1;
    for v = 1:nSlices
        for b = 1:nBaseStations
            opVec(ptr : (ptr + loadDist(v, b) - 1)) = v;
            bsVec(ptr : (ptr + loadDist(v, b) - 1)) = b;
            ptr = ptr + loadDist(v, b);
        end
    end
    bsAssociation{t} = bsVec;
    capacities{t} = capVec;

    % Adjust profile for backward compatibility.
    tmpNetSettings = netSettings;
    tmpNetSettings.users = length(bsAssociation{t});
    if (tmpNetSettings.users == 0)
        continue;
    end
    tmpOpSettings = opSettings;
    tmpOpSettings.s_o = shareVec;
    tmpOpSettings.ops_belongs = opVec;
    tmpOpSettings.shareDist = shareDist;
    weightPerUser = zeros(1, tmpNetSettings.users);
    for v = 1:nSlices
        weightPerUser(tmpOpSettings.ops_belongs == v) = shareVec(v) ...
            ./ sum(tmpOpSettings.ops_belongs == v);
    end
    tmpOpSettings.w_i = weightPerUser;
    perUserMinRateReq = zeros(1, nUsers);
    for v = 1:nSlices
        if (sliceCats(v) == 0)
            perUserMinRateReq(tmpOpSettings.ops_belongs == v) = minRateReq(v);
        end
    end
        
    [r, f, b, nRounds, isViolation] = DIFFPRICE(tmpNetSettings, tmpOpSettings, capVec, ...
        bsVec, minRateReq, 0, sliceCats);

    nRoundsVec(t) = nRounds;
    violationVec(t) = isViolation;
    
    fprintf('finish at time %d\n', t);
end

%% Plot cdf of convergence
datestring = datestr(now, 30);

figure(2)
hold on; grid on
%cdf1 = cdfplot(nRoundsVec(violationVec == 0));
cdf2 = cdfplot(nRoundsVec(violationVec == 1));
%set(cdf1, 'marker', '+');
set(cdf2, 'marker', 'o');
legend('When module < 1', 'When module > 1');
xlabel('number of rounds');
ylabel('CDF');
savefig(sprintf('figs/cdf-convergence-%s.fig', datestring));


%% 
