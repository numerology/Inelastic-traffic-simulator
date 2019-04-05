% Simulation script for user mobility induced by Poisson processes.
% Author: numerology
% 3/15/2019
clc, close all, clear all, gcp;
%% Settings
nSlice = 3;
simulationTime = 5000;
shareVec = 3 * [0.4 0.3 0.3];
perBSLoad = 6;
rhoVec = perBSLoad * [0.4 0.1 0.5;0.4 0.3 0.3;0.4 0.5 0.1]'; % mean load distribution, V x B
shareDist = [0.4 0.1 0.5;0.4 0.3 0.3;0.4 0.5 0.1]';

nBaseStations = size(rhoVec, 2); % Since it's a simpler model, sectors are not mentioned
capacity = 1; % Uniform fixed capacity

threshold = 0.5 * capacity / perBSLoad;

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
    loadDist = 1 * poissrnd(rhoVec);
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
%     [r, f, b] = GPS(tmpNetSettings, tmpOpSettings, capacities{t}, ...
%         bsAssociation{t});
%     ratesStaticGPS{t} = r;
%     [r, f, b] = flexibleSS(tmpNetSettings, tmpOpSettings, capacities{t}, ...
%         bsAssociation{t});
%     ratesSS{t} = r;
%     [r, f, b] = Static_Slicing(tmpNetSettings, tmpOpSettings, capacities{t}, ...
%         bsAssociation{t});
%     ratesStaticSS{t} = r;
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
%flatRateStaticGPS = horzcat(ratesStaticGPS{:});
flatRateGPS = horzcat(ratesGPS{:});
%flatRateStaticSS = horzcat(ratesStaticSS{:});
%flatRateSS = horzcat(ratesSS{:});
flatRateMW = horzcat(ratesMW{:});
flatRateMWBR = horzcat(ratesMWBR{:});
flatRateSCPF = horzcat(ratesSCPF{:});

figure()
hold on
%cdfplot(1./flatRateStaticSS);
%cdfplot(1./flatRateSS);
%cdfplot(1./flatRateStaticGPS);
cdfplot(1./flatRateGPS);
handleMWPA = cdfplot(1./flatRateMW);
cdfplot(1./flatRateMWBR);
handleSCPF = cdfplot(1./flatRateSCPF);
title('CDF of BTD')
xlabel('BTD')
legend('GPS', 'MAXWEIGHT - PA', ...
    'MAXWEIGHT - BR', 'SCPF');
%set(handleSCPF, 'Marker', 'o');
%set(handleMWPA, 'Marker', '+');

figure()
hold on
% cdfplot(flatRateStaticSS);
% cdfplot(flatRateSS);
% cdfplot(flatRateStaticGPS);
cdfplot(flatRateGPS);
handleMWPA = cdfplot(flatRateMW);
cdfplot(flatRateMWBR);
handleSCPF = cdfplot(flatRateSCPF);
title('CDF of rate')
xlabel('rate')
legend('GPS', 'MAXWEIGHT - PA',...
    'MAXWEIGHT - BR', 'SCPF');
%set(handleSCPF, 'Marker', 'o');
%set(handleMWPA, 'Marker', '+');

disp('Overall')
fprintf('mean btd of static SS = %f\n', mean(1 ./ flatRateStaticSS));
fprintf('mean btd of SS = %f\n', mean(1 ./ flatRateSS));
fprintf('mean btd of static GPS = %f\n', mean(1 ./ flatRateStaticGPS));
fprintf('mean btd of GPS = %f\n', mean(1 ./ flatRateGPS));
fprintf('mean btd of SCPF = %f\n', mean(1 ./ flatRateSCPF));
fprintf('mean btd of MAXWEIGHT-PA = %f\n', mean(1 ./ flatRateMW));
fprintf('mean btd of MAXWEIGHT-BR = %f\n', mean(1 ./ flatRateMWBR));
%% Slice specific analysis
ratesStaticGPS1 = filterbyslice(ratesStaticGPS, opBelongs, 1);
ratesGPS1 = filterbyslice(ratesGPS, opBelongs, 1);
ratesStaticSS1 = filterbyslice(ratesStaticSS, opBelongs, 1);
ratesSS1 = filterbyslice(ratesSS, opBelongs, 1);
ratesMW1 = filterbyslice(ratesMW, opBelongs, 1);
ratesMWBR1 = filterbyslice(ratesMWBR, opBelongs, 1);
ratesSCPF1 = filterbyslice(ratesSCPF, opBelongs, 1);

flatRateStaticGPS1 = horzcat(ratesStaticGPS1{:});
flatRateGPS1 = horzcat(ratesGPS1{:});
flatRateStaticSS1 = horzcat(ratesStaticSS1{:});
flatRateSS1 = horzcat(ratesSS1{:});
flatRateMW1 = horzcat(ratesMW1{:});
flatRateMWBR1 = horzcat(ratesMWBR1{:});
flatRateSCPF1 = horzcat(ratesSCPF1{:});

figure()
hold on
cdfplot(1./flatRateStaticSS1);
cdfplot(1./flatRateSS1);
cdfplot(1./flatRateStaticGPS1);
cdfplot(1./flatRateGPS1);
handleMWPA = cdfplot(1./flatRateMW1);
cdfplot(1./flatRateMWBR1);
handleSCPF = cdfplot(1./flatRateSCPF1);
title('CDF of BTD of slice 1')
xlabel('BTD')
legend('Static SS', 'SS', 'Static GPS', 'GPS', 'MAXWEIGHT - PA', ...
    'MAXWEIGHT - BR', 'SCPF');
set(handleSCPF, 'Marker', 'o');
set(handleMWPA, 'Marker', '+');

figure()
hold on
cdfplot(flatRateStaticSS1);
cdfplot(flatRateSS1);
cdfplot(flatRateStaticGPS1);
cdfplot(flatRateGPS1);
handleMWPA = cdfplot(flatRateMW1);
cdfplot(flatRateMWBR1);
handleSCPF = cdfplot(flatRateSCPF1);
title('CDF of rate of slice 1')
xlabel('rate')
legend('Static SS', 'SS', 'Static GPS', 'GPS', 'MAXWEIGHT - PA', ...
    'MAXWEIGHT - BR', 'SCPF');
set(handleSCPF, 'Marker', 'o');
set(handleMWPA, 'Marker', '+');

disp('For a typical user on slice 1')
fprintf('mean btd of static SS = %f\n', mean(1 ./ flatRateStaticSS1));
fprintf('mean btd of SS = %f\n', mean(1 ./ flatRateSS1));
fprintf('mean btd of static GPS = %f\n', mean(1 ./ flatRateStaticGPS1));
fprintf('mean btd of GPS = %f\n', mean(1 ./ flatRateGPS1));
fprintf('mean btd of SCPF = %f\n', mean(1 ./ flatRateSCPF1));
fprintf('mean btd of MAXWEIGHT-PA = %f\n', mean(1 ./ flatRateMW1));
fprintf('mean btd of MAXWEIGHT-BR = %f\n', mean(1 ./ flatRateMWBR1));

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
handleSCPF = cdfplot(utilSCPF);
handleMWPA = cdfplot(utilMW);
cdfplot(utilMWBR);
title('CDF of utility function');
xlabel('utility')
legend('static SS', 'SS', 'static GPS', 'GPS', 'SCPF', 'MAXWEIGHT-PA', ...
    'MAXWEIGHT-BR');
set(handleSCPF, 'Marker', 'o');
set(handleMWPA, 'Marker', '+');
    

