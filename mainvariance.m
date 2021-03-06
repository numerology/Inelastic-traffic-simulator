% mainvariance The script used to evaluate the impact of changing the
% variance of the spatial load distribution of users.
% Author: numerology
% 03/21/2019
clc, close all, clear all, gcp;
%% Settings
nSlice = 3;

simulationTime = 50000;
perBSLoad = 6;

shareVec = [14/9 13/18 13/18];
relativeRhoVec = [perBSLoad * [2/3 1/6 1/6];
                  perBSLoad * [2/3 1/6 1/6];
                  perBSLoad * [2/9 7/18 7/18]]'; % mean load distribution, V x B
shareDist = [2/3 1/6 1/6;
    2/3 1/6 1/6;
    2/9 7/18 7/18]';

nBaseStations = size(relativeRhoVec, 2);
capacity = 1;
threshold = 0.5 * capacity / perBSLoad; % min rate requirement
netSettings = [];
netSettings.bsNS = nBaseStations;
opSettings = [];
opSettings.s_o = shareVec;
opSettings.shareDist = shareDist;

varFactors = 1:5;
btdGainVecSCPF = zeros(1, length(varFactors)); % BTD gain over (flexible) GPS.
btdGainVecMWPA = zeros(1, length(varFactors));
btdGainVecMWBR = zeros(1, length(varFactors));
utilityGainVecSCPF = zeros(1, length(varFactors)); % overall utility gain over (flexible) GPS.
utilityGainVecMWPA = zeros(1, length(varFactors));
utilityGainVecMWBR = zeros(1, length(varFactors));
ratesGPS = cell(length(varFactors), simulationTime); % Save ordinary data for regression.
ratesMW = cell(length(varFactors), simulationTime);
ratesMWBR = cell(length(varFactors), simulationTime);
ratesSCPF = cell(length(varFactors), simulationTime);
opBelongs = cell(length(varFactors), simulationTime);

meanBtdGPS = zeros(1, length(varFactors));
meanBtdMW = zeros(1, length(varFactors));
meanBtdMWBR = zeros(1, length(varFactors));
meanBtdSCPF = zeros(1, length(varFactors));
meanEffRateGPS = zeros(1, length(varFactors));
meanEffRateMW = zeros(1, length(varFactors));
meanEffRateMWBR = zeros(1, length(varFactors));
meanEffRateSCPF = zeros(1, length(varFactors));
meanUtilGPS = zeros(1, length(varFactors));
meanUtilMW = zeros(1, length(varFactors));
meanUtilMWBR = zeros(1, length(varFactors));
meanUtilSCPF = zeros(1, length(varFactors));

%% Run simulations
for i = 1:length(varFactors)
    varFactor = varFactors(i);
    rhoVec = relativeRhoVec / varFactor;
    bsAssociation = cell(1, simulationTime);
    capacities = cell(1, simulationTime);
    ppm = ParforProgMon('simulating: ', simulationTime);
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
        opBelongs{i, t} = opVec;
        capacities{t} = capVec;
        
        % Adjust profile for backward compatibility.
        tmpNetSettings = netSettings;
        tmpNetSettings.users = length(bsAssociation{t});
        if (tmpNetSettings.users == 0)
            ppm.increment();
            continue;
        end
        tmpOpSettings = opSettings;
        tmpOpSettings.ops_belongs = opBelongs{i, t};
        tmpOpSettings.shareDist = shareDist;
        weightPerUser = zeros(1, tmpNetSettings.users);
        for v = 1:nSlice
            weightPerUser(tmpOpSettings.ops_belongs == v) = shareVec(v) ...
                ./ sum(tmpOpSettings.ops_belongs == v);
        end
        tmpOpSettings.w_i = weightPerUser;
        [r, f, b] = flexibleGPS(tmpNetSettings, tmpOpSettings, capacities{t}, ...
            bsAssociation{t});
        ratesGPS{i, t} = r;
        [r, f, b] = SCPF(tmpNetSettings, tmpOpSettings, capacities{t}, ...
            bsAssociation{t});
        ratesSCPF{i, t} = r;
        [r, f, b] = MAXWEIGHT_PA(tmpNetSettings, tmpOpSettings, capacities{t}, ...
            bsAssociation{t});
        ratesMW{i, t} = r;
        [r, f, b] = MAXWEIGHT_BR(tmpNetSettings, tmpOpSettings, capacities{t}, ...
            bsAssociation{t});
        ratesMWBR{i, t} = r;
        if (sum(ratesMWBR{i, t} < 1e-4) > 0)
            ratesMWBR{i, t}(ratesMWBR{i, t} < 1e-4) = nan;
        end
        ppm.increment();
    end
    flatRateGPS = horzcat(ratesGPS{i, :});
    flatRateMW = horzcat(ratesMW{i, :});
    flatRateMWBR = horzcat(ratesMWBR{i, :});
    flatRateSCPF = horzcat(ratesSCPF{i, :});
    meanBtdGPS(i) = mean(1./flatRateGPS);
    meanBtdMW(i) = mean(1./flatRateMW);
    meanBtdMWBR(i) = nanmean(1./flatRateMWBR);
    meanBtdSCPF(i) = mean(1./flatRateSCPF);
  
    btdGainVecSCPF(i) = mean(1./flatRateGPS) / mean(1./flatRateSCPF); 
    btdGainVecMWPA(i) = mean(1./flatRateGPS) / mean(1./flatRateMW);
    btdGainVecMWBR(i) = mean(1./flatRateGPS) ...
        / nanmean(1./flatRateMWBR(flatRateMWBR > 1e-4));
    utilGPS = zeros(1, simulationTime);
    utilSCPF = zeros(1, simulationTime);
    utilMW = zeros(1, simulationTime);
    utilMWBR = zeros(1, simulationTime);
    
    parfor t = 1:simulationTime % stats
        utilGPS(t) = ratetoutil(ratesGPS{i, t}, shareVec, opBelongs{i, t});
        utilSCPF(t) = ratetoutil(ratesSCPF{i, t}, shareVec, opBelongs{i, t});
        utilMW(t) = ratetoutil(ratesMW{i, t}, shareVec, opBelongs{i, t});
        if (sum(ratesMWBR{i, t} < 1e-4) > 0)
            utilMWBR(t) = nan;
        else
            utilMWBR(t) = ratetoutil(ratesMWBR{i, t}, shareVec, opBelongs{i, t});
        end
    end
    meanUtilGPS(i) = nanmean(utilGPS);
    meanUtilMW(i) = nanmean(utilMW);
    meanUtilMWBR(i) = nanmean(utilMWBR);
    meanUtilSCPF(i) = nanmean(utilSCPF);
    utilityGainVecSCPF(i) = nanmean(utilGPS) / nanmean(utilSCPF); 
    utilityGainVecMWPA(i) = nanmean(utilGPS) / nanmean(utilMW);
    utilityGainVecMWBR(i) = nanmean(utilGPS) / nanmean(utilMWBR);
    
    flatRateGPS(flatRateGPS < threshold) = 0;
    flatRateMW(flatRateMW < threshold) = 0;
    flatRateMWBR(flatRateMWBR < threshold) = 0;
    flatRateSCPF(flatRateSCPF < threshold) = 0;
    meanEffRateGPS(i) = nanmean(flatRateGPS);
    meanEffRateMW(i) = nanmean(flatRateMW);
    meanEffRateMWBR(i) = nanmean(flatRateMWBR);
    meanEffRateSCPF(i) = nanmean(flatRateSCPF);
end

figure(7);
hold on
plot(varFactors, meanEffRateGPS, 'gd-');
plot(varFactors, meanEffRateSCPF, 'b+-');
plot(varFactors, meanEffRateMW, 'ro-');
plot(varFactors, meanEffRateMWBR, 'kx-');
title('Average rate above threshold vs. variance factor');
legend('GPS', 'SCPF', 'MAXWEIGHT-practical approach', 'MAXWEIGHT-best response');
saveas(gcf, 'effrate-vs-var.fig');

figure(8);
hold on
plot(varFactors, meanEffRateSCPF./meanEffRateGPS, 'b+-');
plot(varFactors, meanEffRateMW./meanEffRateGPS, 'ro-');
plot(varFactors, meanEffRateMWBR./meanEffRateGPS, 'kx-');
title('Average rate above threshold gain over GPS vs. variance factor');
legend('SCPF', 'MAXWEIGHT-practical approach', 'MAXWEIGHT-best response');
saveas(gcf, 'effrate-gain-vs-var.fig');

figure(1);
hold on
plot(varFactors, btdGainVecSCPF, 'b+-');
plot(varFactors, btdGainVecMWPA, 'ro-');
plot(varFactors, btdGainVecMWBR, 'kx-');
title('BTD gain over GPS vs. variance factor');
legend('SCPF', 'MAXWEIGHT-practical approach', 'MAXWEIGHT-best response');
saveas(gcf, 'btd-gain-vs-var.fig');

figure(2);
hold on
plot(varFactors, meanUtilSCPF - meanUtilGPS, 'b+-');
plot(varFactors, meanUtilMW - meanUtilGPS, 'ro-');
plot(varFactors, meanUtilMWBR - meanUtilGPS, 'kx-');
title('Utility gain over GPS vs. variance factor');
legend('SCPF', 'MAXWEIGHT-practical approach', 'MAXWEIGHT-best response');
saveas(gcf, 'util-gain-vs-var.fig');

% Absolute rate
figure(3);
hold on
plot(varFactors, meanBtdGPS, 'gd-');
plot(varFactors, meanBtdSCPF, 'b+-');
plot(varFactors, meanBtdMW, 'ro-');
plot(varFactors, meanBtdMWBR, 'kx-');
title('Average btd vs. variance factor');
legend('GPS', 'SCPF', 'MAXWEIGHT-practical approach', 'MAXWEIGHT-best response');
saveas(gcf, 'btd-vs-var.fig');

% Absolute utility
figure(4);
hold on
plot(varFactors, meanUtilGPS, 'gd-');
plot(varFactors, meanUtilSCPF, 'b+-');
plot(varFactors, meanUtilMW, 'ro-');
plot(varFactors, meanUtilMWBR, 'kx-');
title('Average utility vs. variance factor');
legend('GPS', 'SCPF', 'MAXWEIGHT-practical approach', 'MAXWEIGHT-best response');
saveas(gcf, 'util-vs-var.fig');

%% Get some idea on slice idx
idx = 1;
btdGainVecSCPF1 = zeros(1, length(varFactors)); % BTD gain over (flexible) GPS.
btdGainVecMWPA1 = zeros(1, length(varFactors));
btdGainVecMWBR1 = zeros(1, length(varFactors));
utilityGainVecSCPF1 = zeros(1, length(varFactors)); % overall utility gain over (flexible) GPS.
utilityGainVecMWPA1 = zeros(1, length(varFactors));
utilityGainVecMWBR1 = zeros(1, length(varFactors));
meanBtdGPS1 = zeros(1, length(varFactors));
meanBtdMW1 = zeros(1, length(varFactors));
meanBtdMWBR1 = zeros(1, length(varFactors));
meanBtdSCPF1 = zeros(1, length(varFactors));
meanUtilGPS1 = zeros(1, length(varFactors));
meanUtilMW1 = zeros(1, length(varFactors));
meanUtilMWBR1 = zeros(1, length(varFactors));
meanUtilSCPF1 = zeros(1, length(varFactors));

for i = 1:length(varFactors)
    sliceIdx = (horzcat(opBelongs{i, :}) == idx);
    flatRateGPS1 = horzcat(ratesGPS{i, :});
    flatRateGPS1 = flatRateGPS1(sliceIdx);
    flatRateMW1 = horzcat(ratesMW{i, :});
    flatRateMW1 = flatRateMW1(sliceIdx);
    flatRateMWBR1 = horzcat(ratesMWBR{i, :});
    flatRateMWBR1 = flatRateMWBR1(sliceIdx);
    flatRateMWBR1 = flatRateMWBR1(flatRateMWBR1 > 1e-4);
    flatRateSCPF1 = horzcat(ratesSCPF{i, :});
    flatRateSCPF1 = flatRateSCPF1(sliceIdx);
    meanBtdGPS1(i) = mean(1./flatRateGPS1);
    meanBtdMW1(i) = mean(1./flatRateMW1);
    meanBtdMWBR1(i) = nanmean(1./flatRateMWBR1);
    meanBtdSCPF1(i) = mean(1./flatRateSCPF1);
    
    btdGainVecSCPF1(i) = mean(1./flatRateGPS1) / mean(1./flatRateSCPF1); 
    btdGainVecMWPA1(i) = mean(1./flatRateGPS1) / mean(1./flatRateMW1);
    btdGainVecMWBR1(i) = mean(1./flatRateGPS1) ...
        / nanmean(1./flatRateMWBR1);
    
    % (TODO:) Make utility computation per slice.
%     utilGPS1 = zeros(1, simulationTime);
%     utilSCPF1 = zeros(1, simulationTime);
%     utilMW1 = zeros(1, simulationTime);
%     utilMWBR1 = zeros(1, simulationTime);
% 
%     parfor t = 1:simulationTime
%         utilGPS(t) = ratetoutil(ratesGPS{i, t}, shareVec, opBelongs{i, t});
%         utilSCPF(t) = ratetoutil(ratesSCPF{i, t}, shareVec, opBelongs{i, t});
%         utilMW(t) = ratetoutil(ratesMW{i, t}, shareVec, opBelongs{i, t});
%         utilMWBR(t) = ratetoutil(ratesMWBR{i, t}, shareVec, opBelongs{i, t});
%     end
%     meanUtilGPS(i) = nanmean(utilGPS);
%     meanUtilMW(i) = nanmean(utilMW);
%     meanUtilMWBR(i) = nanmean(utilMWBR);
%     meanUtilSCPF(i) = nanmean(utilSCPF);
%     utilityGainVecSCPF(i) = nanmean(utilGPS) / nanmean(utilSCPF); 
%     utilityGainVecMWPA(i) = nanmean(utilGPS) / nanmean(utilMW);
%     utilityGainVecMWBR(i) = nanmean(utilGPS) / nanmean(utilMWBR);
end

figure(5)
hold on
plot(varFactors, meanBtdGPS1, 'gd-');
plot(varFactors, meanBtdSCPF1, 'b+-');
plot(varFactors, meanBtdMW1, 'ro-');
plot(varFactors, meanBtdMWBR1, 'kx-');
title('Average btd vs. variance factor of slice 1');
legend('GPS', 'SCPF', 'MAXWEIGHT-practical approach', 'MAXWEIGHT-best response');

figure(6)
hold on
plot(varFactors, btdGainVecSCPF1, 'b+-');
plot(varFactors, btdGainVecMWPA1, 'ro-');
plot(varFactors, btdGainVecMWBR1, 'kx-');
title('BTD gain over GPS vs. variance factor on slice 1');
legend('SCPF', 'MAXWEIGHT-practical approach', 'MAXWEIGHT-best response');