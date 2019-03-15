% Simulation script for user mobility induced by Poisson processes.
% Author: numerology
% 3/15/2019
clc, close all, clear all, gcp;
%% Settings
nSlice = 3;
simulationTime = 1000;
shareVec = 1/3 * ones(1, 3);
rhoVec = 100 * [0.5 0.25 0.25;0.25 0.5 0.25;0.25 0.25 0.5]'; % mean load distribution, V x B
shareDist = [0.5 0.25 0.25;0.25 0.5 0.25;0.25 0.25 0.5]';

nBaseStations = size(rhoVec, 2); % Since it's a simpler model, sectors are not mentioned
capacity = 1; % Uniform fixed capacity

netSettings = [];
netSettings.bsNS = nBaseStations;
opSettings = [];
opSettings.s_o = shareVec;
opSettings.shareDist = shareDist;

%% Generate the user association
bsAssociation = cell(1, simulationTime);
opBelongs = cell(1, simulationTime);
capacities = cell(1, simulationTime);
ppm = ParforProgMon('Generating network profile: ', simulationTime);
parfor t = 1:simulationTime
    loadDist = poissrnd(rhoVec);
    nUsers = sum(sum(loadDist));
    bsVec = zeros(1, nUsers);
    opVec = zeros(1, nUsers);
    capVec = capacity * ones(1, nUsers);
    ptr = 1;
    for v = 1:nSlice
        for b = 1:nBaseStations
            opVec(ptr : (ptr + loadDist(v, b) - 1)) = v;
            bsVec(ptr : (ptr + loadDist(v, b) - 1)) = b;
            ptr = ptr + loadDist(v, b);
        end
    end
    bsAssociation{t} = bsVec;
    opBelongs{t} = opVec;
    capacities{t} = capVec;
    ppm.increment();
end

%% Run different sharing criteria
ratesGPS = cell(1, simulationTime);
ratesMW = cell(1, simulationTime);
ratesSCPF = cell(1, simulationTime);
ppm = ParforProgMon('Simulating resource sharing : ', simulationTime);
parfor t = 1:simulationTime
    % Adjust profile for backward compatibility.
    tmpNetSettings = netSettings;
    tmpNetSettings.users = length(bsAssociation{t});
    if (tmpNetSettings.users == 0)
        ppm.increment();
        continue;
    end
    tmpOpSettings = opSettings;
    tmpOpSettings.ops_belongs = opBelongs{t};
    tmpOpSettings.shareDist = shareDist;
    weightPerUser = zeros(1, tmpNetSettings.users);
    for v = 1:nSlice
        weightPerUser(tmpOpSettings.ops_belongs == v) = shareVec(v) ./ sum(...
            tmpOpSettings.ops_belongs == v);
    end
    tmpOpSettings.w_i = weightPerUser;
    [r, f, b] = flexibleGPS(tmpNetSettings, tmpOpSettings, capacities{t}, ...
        bsAssociation{t});
    ratesGPS{t} = r;
    [r, f, b] = SCPF(tmpNetSettings, tmpOpSettings, capacities{t}, ...
        bsAssociation{t});
    ratesSCPF{t} = r;
    [r, f, b] = MAXWEIGHT(tmpNetSettings, tmpOpSettings, capacities{t}, ...
        bsAssociation{t});
    ratesMW{t} = r;
    ppm.increment();
end

%% Some analysis
flatRateGPS = horzcat(ratesGPS{:});
flatRateMW = horzcat(ratesMW{:});
flatRateSCPF = horzcat(ratesSCPF{:});

figure()
hold on
cdfplot(1./flatRateGPS);
cdfplot(1./flatRateMW);
cdfplot(1./flatRateSCPF);
title('CDF of BTD')
xlabel('BTD')
legend('GPS', 'MAXWEIGHT', 'SCPF');

disp('Overall')
fprintf('mean btd of GPS = %f\n', mean(1 ./ flatRateGPS));
fprintf('mean btd of SCPF = %f\n', mean(1 ./ flatRateSCPF));
fprintf('mean btd of MAXWEIGHT = %f\n', mean(1 ./ flatRateMW));


