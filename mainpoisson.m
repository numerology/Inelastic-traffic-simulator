% Simulation script for user mobility induced by Poisson processes.
% Author: numerology
% 3/15/2019
clc, close all, clear all, gcp;
%% Settings
nSlice = 3;
simulationTime = 200;
shareVec = 1/3 * ones(1, 3);
rhoVec = 10 * [0.8 0.1 0.1;0.1 0.8 0.1;0.1 0.1 0.8]'; % mean load distribution, V x B
shareDist = [0.8 0.1 0.1;0.1 0.8 0.1;0.1 0.1 0.8]';

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
ratesStaticGPS = cell(1, simulationTime);
ratesSS = cell(1, simulationTime);
ratesStaticSS = cell(1, simulationTime);
ratesMW = cell(1, simulationTime);
ratesMWBR = cell(1, simulationTime);
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
    % Resource sharing
    [r, f, b] = flexibleGPS(tmpNetSettings, tmpOpSettings, capacities{t}, ...
        bsAssociation{t});
    ratesGPS{t} = r;
    [r, f, b] = GPS(tmpNetSettings, tmpOpSettings, capacities{t}, ...
        bsAssociation{t});
    ratesStaticGPS{t} = r;
    [r, f, b] = flexibleSS(tmpNetSettings, tmpOpSettings, capacities{t}, ...
        bsAssociation{t});
    ratesSS{t} = r;
    [r, f, b] = Static_Slicing(tmpNetSettings, tmpOpSettings, capacities{t}, ...
        bsAssociation{t});
    ratesStaticSS{t} = r;
    [r, f, b] = SCPF(tmpNetSettings, tmpOpSettings, capacities{t}, ...
        bsAssociation{t});
    ratesSCPF{t} = r;
    [r, f, b] = MAXWEIGHT_PA(tmpNetSettings, tmpOpSettings, capacities{t}, ...
        bsAssociation{t});
    ratesMW{t} = r;
    [r, f, b] = MAXWEIGHT_BR(tmpNetSettings, tmpOpSettings, capacities{t}, ...
        bsAssociation{t});
    ratesMWBR{t} = r;
    ppm.increment();
end

%% Some analysis
flatRateStaticGPS = horzcat(ratesStaticGPS{:});
flatRateGPS = horzcat(ratesGPS{:});
flatRateStaticSS = horzcat(ratesStaticSS{:});
flatRateSS = horzcat(ratesSS{:});
flatRateMW = horzcat(ratesMW{:});
flatRateMWBR = horzcat(ratesMWBR{:});
flatRateSCPF = horzcat(ratesSCPF{:});

figure()
hold on
cdfplot(1./flatRateStaticSS);
cdfplot(1./flatRateSS);
cdfplot(1./flatRateStaticGPS);
cdfplot(1./flatRateGPS);
cdfplot(1./flatRateMW);
cdfplot(1./flatRateMWBR);
cdfplot(1./flatRateSCPF);
title('CDF of BTD')
xlabel('BTD')
legend('Static SS', 'SS', 'Static GPS', 'GPS', 'MAXWEIGHT - PA', ...
    'MAXWEIGHT - BR', 'SCPF');

disp('Overall')
fprintf('mean btd of static SS = %f\n', mean(1 ./ flatRateStaticSS));
fprintf('mean btd of SS = %f\n', mean(1 ./ flatRateSS));
fprintf('mean btd of static GPS = %f\n', mean(1 ./ flatRateStaticGPS));
fprintf('mean btd of GPS = %f\n', mean(1 ./ flatRateGPS));
fprintf('mean btd of SCPF = %f\n', mean(1 ./ flatRateSCPF));
fprintf('mean btd of MAXWEIGHT-PA = %f\n', mean(1 ./ flatRateMW));
fprintf('mean btd of MAXWEIGHT-BR = %f\n', mean(1 ./ flatRateMWBR));

%% Exam the utility function

utilStaticSS = zeros(1, simulationTime);
utilSS = zeros(1, simulationTime);
utilStaticGPS = zeros(1, simulationTime);
utilGPS = zeros(1, simulationTime);
utilSCPF = zeros(1, simulationTime);
utilMW = zeros(1, simulationTime);
utilMWBR = zeros(1, simulationTime);

for t = 1:simulationTime
    utilStaticSS(t) = ratetoutil(ratesStaticSS{t}, shareVec, opBelongs{t});
    utilSS(t) = ratetoutil(ratesSS{t}, shareVec, opBelongs{t});
    utilStaticGPS(t) = ratetoutil(ratesStaticGPS{t}, shareVec, opBelongs{t});
    utilGPS(t) = ratetoutil(ratesGPS{t}, shareVec, opBelongs{t});
    utilSCPF(t) = ratetoutil(ratesSCPF{t}, shareVec, opBelongs{t});
    utilMW(t) = ratetoutil(ratesMW{t}, shareVec, opBelongs{t});
    utilMWBR(t) = ratetoutil(ratesMWBR{t}, shareVec, opBelongs{t});
end

figure()
hold on
cdfplot(utilStaticSS);
cdfplot(utilSS);
cdfplot(utilStaticGPS);
cdfplot(utilGPS);
cdfplot(utilSCPF);
cdfplot(utilMW);
cdfplot(utilMWBR);
title('CDF of utility function');
xlabel('utility')
legend('static SS', 'SS', 'static GPS', 'GPS', 'SCPF', 'MAXWEIGHT-PA', ...
    'MAXWEIGHT-BR');
    

