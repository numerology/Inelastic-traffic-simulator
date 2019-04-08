% no gui version of mainminreq.m
clc, close all, clear all, gcp;
nSlice = 3;

simulationTime = 10000;
perBSLoad = 6;
shareVec = [14/9 13/18 13/18];
relativeRhoVec = [perBSLoad * [2/3 1/6 1/6];
                  perBSLoad * [2/3 1/6 1/6];
                  perBSLoad * [2/9 7/18 7/18]]'; % mean load distribution, V x B

nBaseStations = size(relativeRhoVec, 2);
capacity = 1;
minRateReq = 0.1 * capacity / perBSLoad * ones(1, nSlice); % min rate requirement
minSharePerBS = 0.1;
outageTol = 0.2;
netSettings = [];
netSettings.bsNS = nBaseStations;
opSettings = [];
opSettings.s_o = shareVec;

varFactors = 1:5;
btdGainVecSCPF = zeros(1, length(varFactors)); % BTD gain over (flexible) GPS.
btdGainVecDP = zeros(1, length(varFactors));
btdGainVecMWBR = zeros(1, length(varFactors));
utilityGainVecSCPF = zeros(1, length(varFactors)); % overall utility gain over (flexible) GPS.
utilityGainVecDP = zeros(1, length(varFactors));
utilityGainVecMWBR = zeros(1, length(varFactors));
ratesGPS = cell(length(varFactors), simulationTime); % Save ordinary data for regression.
ratesDP = cell(length(varFactors), simulationTime);
ratesMWBR = cell(length(varFactors), simulationTime);
ratesSCPF = cell(length(varFactors), simulationTime);
opBelongs = cell(length(varFactors), simulationTime);

meanBtdGPS = zeros(1, length(varFactors));
meanBtdDP = zeros(1, length(varFactors));
meanBtdMWBR = zeros(1, length(varFactors));
meanBtdSCPF = zeros(1, length(varFactors));
meanEffRateGPS = zeros(1, length(varFactors));
meanEffRateDP = zeros(1, length(varFactors));
meanEffRateMWBR = zeros(1, length(varFactors));
meanEffRateSCPF = zeros(1, length(varFactors));
meanUtilGPS = zeros(1, length(varFactors));
meanUtilDP = zeros(1, length(varFactors));
meanUtilMWBR = zeros(1, length(varFactors));
meanUtilSCPF = zeros(1, length(varFactors));

%% Run simulations
for i = 1:length(varFactors)
    varFactor = varFactors(i);
    rhoVec = relativeRhoVec / varFactor;
    bsAssociation = cell(1, simulationTime);
    capacities = cell(1, simulationTime);
    shareDist = sharedimension(minRateReq, rhoVec, shareVec, outageTol, minSharePerBS);
    
    parfor t = 1:simulationTime
        loadDist = ceil(varFactor * poissrnd(rhoVec));
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
        [r, f, b] = DIFFPRICE(tmpNetSettings, tmpOpSettings, capacities{t}, ...
            bsAssociation{t}, minRateReq);
        ratesDP{i, t} = r;
        [r, f, b] = MAXWEIGHT_BR(tmpNetSettings, tmpOpSettings, capacities{t}, ...
            bsAssociation{t});
        ratesMWBR{i, t} = r;
        if (sum(ratesMWBR{i, t} < 1e-4) > 0)
            ratesMWBR{i, t}(ratesMWBR{i, t} < 1e-4) = nan;
        end
    end
    flatRateGPS = horzcat(ratesGPS{i, :});
    flatRateDP = horzcat(ratesDP{i, :});
    flatRateMWBR = horzcat(ratesMWBR{i, :});
    flatRateSCPF = horzcat(ratesSCPF{i, :});
    meanBtdGPS(i) = mean(1./flatRateGPS);
    meanBtdDP(i) = mean(1./flatRateDP);
    meanBtdMWBR(i) = nanmean(1./flatRateMWBR);
    meanBtdSCPF(i) = mean(1./flatRateSCPF);
  
    btdGainVecSCPF(i) = mean(1./flatRateGPS) / mean(1./flatRateSCPF); 
    btdGainVecDP(i) = mean(1./flatRateGPS) / mean(1./flatRateDP);
    btdGainVecMWBR(i) = mean(1./flatRateGPS) ...
        / nanmean(1./flatRateMWBR(flatRateMWBR > 1e-4));
    utilGPS = zeros(1, simulationTime);
    utilSCPF = zeros(1, simulationTime);
    utilDP = zeros(1, simulationTime);
    utilMWBR = zeros(1, simulationTime);
    
    parfor t = 1:simulationTime % stats
        utilGPS(t) = ratetoutil(ratesGPS{i, t}, shareVec, opBelongs{i, t});
        utilSCPF(t) = ratetoutil(ratesSCPF{i, t}, shareVec, opBelongs{i, t});
        utilDP(t) = ratetoutil(ratesDP{i, t}, shareVec, opBelongs{i, t});
        if (sum(ratesMWBR{i, t} < 1e-4) > 0)
            utilMWBR(t) = nan;
        else
            utilMWBR(t) = ratetoutil(ratesMWBR{i, t}, shareVec, opBelongs{i, t});
        end
    end
    meanUtilGPS(i) = nanmean(utilGPS);
    meanUtilDP(i) = nanmean(utilDP);
    meanUtilMWBR(i) = nanmean(utilMWBR);
    meanUtilSCPF(i) = nanmean(utilSCPF);
    utilityGainVecSCPF(i) = nanmean(utilGPS) / nanmean(utilSCPF); 
    utilityGainVecDP(i) = nanmean(utilGPS) / nanmean(utilDP);
    utilityGainVecMWBR(i) = nanmean(utilGPS) / nanmean(utilMWBR);
    
end

% Plot results
datestring = datestr(now, 30);

figure(1);
hold on
plot(varFactors, btdGainVecSCPF, 'b+-');
plot(varFactors, btdGainVecDP, 'ro-');
plot(varFactors, btdGainVecMWBR, 'kx-');
title('BTD gain over GPS vs. variance factor');
legend('SCPF', 'DIFFPRICE', 'MAXWEIGHT-best response');
savefig(sprintf('btd-gain-vs-var-%s.fig', datestring));

figure(2);
hold on
plot(varFactors, meanUtilSCPF - meanUtilGPS, 'b+-');
plot(varFactors, meanUtilDP - meanUtilGPS, 'ro-');
plot(varFactors, meanUtilMWBR - meanUtilGPS, 'kx-');
title('Utility gain over GPS vs. variance factor');
legend('SCPF', 'DIFFPRICE', 'MAXWEIGHT-best response');
savefig(sprintf('util-gain-vs-var-%s.fig', datestring));

% Get some idea on slice idx
idx = 2;
btdGainVecSCPF1 = zeros(1, length(varFactors)); % BTD gain over (flexible) GPS.
btdGainVecDP1 = zeros(1, length(varFactors));
btdGainVecMWBR1 = zeros(1, length(varFactors));
meanBtdGPS1 = zeros(1, length(varFactors));
meanBtdDP1 = zeros(1, length(varFactors));
meanBtdMWBR1 = zeros(1, length(varFactors));
meanBtdSCPF1 = zeros(1, length(varFactors));

for i = 1:length(varFactors)
    sliceIdx = (horzcat(opBelongs{i, :}) == idx);
    flatRateGPS1 = horzcat(ratesGPS{i, :});
    flatRateGPS1 = flatRateGPS1(sliceIdx);
    flatRateDP1 = horzcat(ratesDP{i, :});
    flatRateDP1 = flatRateDP1(sliceIdx);
    flatRateMWBR1 = horzcat(ratesMWBR{i, :});
    flatRateMWBR1 = flatRateMWBR1(sliceIdx);
    flatRateMWBR1 = flatRateMWBR1(flatRateMWBR1 > 1e-4);
    flatRateSCPF1 = horzcat(ratesSCPF{i, :});
    flatRateSCPF1 = flatRateSCPF1(sliceIdx);
    meanBtdGPS1(i) = mean(1./flatRateGPS1);
    meanBtdDP1(i) = mean(1./flatRateDP1);
    meanBtdMWBR1(i) = nanmean(1./flatRateMWBR1);
    meanBtdSCPF1(i) = mean(1./flatRateSCPF1);
    
    btdGainVecSCPF1(i) = mean(1./flatRateGPS1) / mean(1./flatRateSCPF1); 
    btdGainVecDP1(i) = mean(1./flatRateGPS1) / mean(1./flatRateDP1);
    btdGainVecMWBR1(i) = mean(1./flatRateGPS1) ...
        / nanmean(1./flatRateMWBR1);
end

figure(3)
hold on
plot(varFactors, meanBtdGPS1, 'gd-');
plot(varFactors, meanBtdSCPF1, 'b+-');
plot(varFactors, meanBtdDP1, 'ro-');
plot(varFactors, meanBtdMWBR1, 'kx-');
title('Average btd vs. variance factor of slice 1');
legend('GPS', 'SCPF', 'MAXWEIGHT-practical approach', 'MAXWEIGHT-best response');
savefig(sprintf('btd-vs-var-slice1-%s.fig', datestring));

figure(4)
hold on
plot(varFactors, btdGainVecSCPF1, 'b+-');
plot(varFactors, btdGainVecDP1, 'ro-');
plot(varFactors, btdGainVecMWBR1, 'kx-');
title('BTD gain over GPS vs. variance factor on slice 1');
legend('SCPF', 'MAXWEIGHT-practical approach', 'MAXWEIGHT-best response');
savefig(sprintf('btd-gain-vs-var-slice1-%s.fig', datestring));

idx = 1;
btdGainVecSCPF1 = zeros(1, length(varFactors)); % BTD gain over (flexible) GPS.
btdGainVecDP1 = zeros(1, length(varFactors));
btdGainVecMWBR1 = zeros(1, length(varFactors));
utilityGainVecSCPF1 = zeros(1, length(varFactors)); % overall utility gain over (flexible) GPS.
utilityGainVecDP1 = zeros(1, length(varFactors));
utilityGainVecMWBR1 = zeros(1, length(varFactors));
meanBtdGPS1 = zeros(1, length(varFactors));
meanBtdDP1 = zeros(1, length(varFactors));
meanBtdMWBR1 = zeros(1, length(varFactors));
meanBtdSCPF1 = zeros(1, length(varFactors));
meanUtilGPS1 = zeros(1, length(varFactors));
meanUtilDP1 = zeros(1, length(varFactors));
meanUtilMWBR1 = zeros(1, length(varFactors));
meanUtilSCPF1 = zeros(1, length(varFactors));

for i = 1:length(varFactors)
    sliceIdx = (horzcat(opBelongs{i, :}) == idx);
    flatRateGPS1 = horzcat(ratesGPS{i, :});
    flatRateGPS1 = flatRateGPS1(sliceIdx);
    flatRateDP1 = horzcat(ratesDP{i, :});
    flatRateDP1 = flatRateDP1(sliceIdx);
    flatRateMWBR1 = horzcat(ratesMWBR{i, :});
    flatRateMWBR1 = flatRateMWBR1(sliceIdx);
    flatRateMWBR1 = flatRateMWBR1(flatRateMWBR1 > 1e-4);
    flatRateSCPF1 = horzcat(ratesSCPF{i, :});
    flatRateSCPF1 = flatRateSCPF1(sliceIdx);
    meanBtdGPS1(i) = mean(1./flatRateGPS1);
    meanBtdDP1(i) = mean(1./flatRateDP1);
    meanBtdMWBR1(i) = nanmean(1./flatRateMWBR1);
    meanBtdSCPF1(i) = mean(1./flatRateSCPF1);
    
    btdGainVecSCPF1(i) = mean(1./flatRateGPS1) / mean(1./flatRateSCPF1); 
    btdGainVecDP1(i) = mean(1./flatRateGPS1) / mean(1./flatRateDP1);
    btdGainVecMWBR1(i) = mean(1./flatRateGPS1) ...
        / nanmean(1./flatRateMWBR1);
end

figure(5)
hold on
plot(varFactors, meanBtdGPS1, 'gd-');
plot(varFactors, meanBtdSCPF1, 'b+-');
plot(varFactors, meanBtdDP1, 'ro-');
plot(varFactors, meanBtdMWBR1, 'kx-');
title('Average btd vs. variance factor of slice 1');
legend('GPS', 'SCPF', 'MAXWEIGHT-practical approach', 'MAXWEIGHT-best response');
savefig(sprintf('btd-vs-var-slice1-%s.fig', datestring));

figure(6)
hold on
plot(varFactors, btdGainVecSCPF1, 'b+-');
plot(varFactors, btdGainVecDP1, 'ro-');
plot(varFactors, btdGainVecMWBR1, 'kx-');
title('BTD gain over GPS vs. variance factor on slice 1');
legend('SCPF', 'MAXWEIGHT-practical approach', 'MAXWEIGHT-best response');
savefig(sprintf('btd-gain-vs-var-slice1-%s.fig', datestring));




