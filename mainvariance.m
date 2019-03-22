% mainvariance The script used to evaluate the impact of changing the
% variance of the spatial load distribution of users.
% Author: numerology
% 03/21/2019
clc, close all, clear all, gcp;
%% Settings
nSlice = 3;
simulationTime = 20000;
shareVec = 3 * [0.4 0.3 0.3];
relativeRhoVec = 3 * [0.4 0.1 0.5;
    0.4 0.3 0.3;
    0.4 0.5 0.1]'; % mean load distribution, V x B
shareDist = [0.4 0.1 0.5;
    0.4 0.3 0.3;
    0.4 0.5 0.1]';
nBaseStations = size(relativeRhoVec, 2);
capacity = 1;

netSettings = [];
netSettings.bsNS = nBaseStations;
opSettings = [];
opSettings.s_o = shareVec;
opSettings.shareDist = shareDist;

varFactors = 1:5;
btdGainVecSCPF = zeros(1, 5); % BTD gain over (flexible) GPS.
btdGainVecMWPA = zeros(1, 5);
btdGainVecMWBR = zeros(1, 5);
utilityGainVecSCPF = zeros(1, 5); % overall utility gain over (flexible) GPS.
utilityGainVecMWPA = zeros(1, 5);
utilityGainVecMWBR = zeros(1, 5);
%% Run simulations
for i = 1:length(varFactors)
    varFactor = varFactors(i);
    rhoVec = relativeRhoVec / varFactor;
    bsAssociation = cell(1, simulationTime);
    opBelongs = cell(1, simulationTime);
    capacities = cell(1, simulationTime);
    ppm = ParforProgMon('Generating network profile: ', simulationTime);
    ratesGPS = cell(1, simulationTime);
    ratesMW = cell(1, simulationTime);
    ratesMWBR = cell(1, simulationTime);
    ratesSCPF = cell(1, simulationTime);
    
    parfor t = 1:simulationTime
        loadDist = varFactor * poissrnd(rhoVec);
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
            weightPerUser(tmpOpSettings.ops_belongs == v) = shareVec(v) ...
                ./ sum(tmpOpSettings.ops_belongs == v);
        end
        tmpOpSettings.w_i = weightPerUser;
        [r, f, b] = flexibleGPS(tmpNetSettings, tmpOpSettings, capacities{t}, ...
            bsAssociation{t});
        ratesGPS{t} = r;
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
    flatRateGPS = horzcat(ratesGPS{:});
    flatRateMW = horzcat(ratesMW{:});
    flatRateMWBR = horzcat(ratesMWBR{:});
    flatRateSCPF = horzcat(ratesSCPF{:});
    btdGainVecSCPF(i) = mean(1./flatRateGPS) / mean(1./flatRateSCPF); 
    btdGainVecMWPA(i) = mean(1./flatRateGPS) / mean(1./flatRateMW);
    btdGainVecMWBR(i) = mean(1./flatRateGPS) ...
        / mean(1./flatRateMWBR(flatRateMWBR > 1e-5));
    utilGPS = zeros(1, simulationTime);
    utilSCPF = zeros(1, simulationTime);
    utilMW = zeros(1, simulationTime);
    utilMWBR = zeros(1, simulationTime);

    parfor t = 1:simulationTime
        utilGPS(t) = ratetoutil(ratesGPS{t}, shareVec, opBelongs{t});
        utilSCPF(t) = ratetoutil(ratesSCPF{t}, shareVec, opBelongs{t});
        utilMW(t) = ratetoutil(ratesMW{t}, shareVec, opBelongs{t});
        utilMWBR(t) = ratetoutil(ratesMWBR{t}, shareVec, opBelongs{t});
    end
    
    utilityGainVecSCPF(i) = nanmean(utilGPS) / nanmean(utilSCPF); 
    utilityGainVecMWPA(i) = nanmean(utilGPS) / nanmean(utilMW);
    utilityGainVecMWBR(i) = nanmean(utilGPS) / nanmean(utilMWBR);
end
%% Plot
figure(1);
hold on
plot(varFactors, btdGainVecSCPF, 'b+-');
plot(varFactors, btdGainVecMWPA, 'ro-');
plot(varFactors, btdGainVecMWBR, 'kx-');
title('BTD gain over GPS vs. variance factor');
legend('SCPF', 'MAXWEIGHT-practical approach', 'MAXWEIGHT-best response');

figure(2);
hold on
plot(varFactors, utilityGainVecSCPF, 'b+-');
plot(varFactors, utilityGainVecMWPA, 'ro-');
plot(varFactors, utilityGainVecMWBR, 'kx-');
title('Utility gain over GPS vs. variance factor');
legend('SCPF', 'MAXWEIGHT-practical approach', 'MAXWEIGHT-best response');
