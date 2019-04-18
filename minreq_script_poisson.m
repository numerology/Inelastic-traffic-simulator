% no gui version of mainminreq.m
clc, close all, clear all;
parpool('local', 40);
warning('off','all');
nSlice = 3;

simulationTime = 1000;
perBSLoad = 6;
% shareVec = [1 1 1];
% relativeRhoVec = [perBSLoad * [1/3 1/3 1/3];
%                   perBSLoad * [1/3 1/3 1/3];
%                   perBSLoad * [1/3 1/3 1/3]]'; % mean load distribution, V x B
shareVec = [1 1 1];
relativeRhoVec = [perBSLoad * [2/3 1/6 1/6];
                  perBSLoad * [1/6 2/3 1/6];
                  perBSLoad * [1/6 1/6 2/3]]'; % mean load distribution, V x B             
% shareVec = [14/9 13/18 13/18];
% relativeRhoVec = [perBSLoad * [2/3 1/6 1/6];
%                   perBSLoad * [2/3 1/6 1/6];
%                   3 * perBSLoad * [2/9 7/18 7/18]]'; % mean load distribution, V x B

nBaseStations = size(relativeRhoVec, 2);
capacity = 1;
minSharePerBS = 0.05;
outageTol = 0.2;
netSettings = [];
netSettings.bsNS = nBaseStations;
opSettings = [];
opSettings.s_o = shareVec;

varFactors = 1:0.5:3;
btdGainVecSCPF = zeros(1, length(varFactors)); % BTD gain over (flexible) GPS.
btdGainVecDP = zeros(1, length(varFactors));
btdGainVecDPWF = zeros(1, length(varFactors));
btdGainVecDPoptimal = zeros(1, length(varFactors));
btdGainVecMWBR = zeros(1, length(varFactors));
utilityGainVecSCPF = zeros(1, length(varFactors)); % overall utility gain over (flexible) GPS.
utilityGainVecDP = zeros(1, length(varFactors));
utilityGainVecDPWF = zeros(1, length(varFactors));
utilityGainVecDPoptimal = zeros(1, length(varFactors));
utilityGainVecMWBR = zeros(1, length(varFactors));
ratesGPS = cell(length(varFactors), simulationTime); % Save ordinary data for regression.
ratesDP = cell(length(varFactors), simulationTime);
ratesDPWF = cell(length(varFactors), simulationTime);
ratesDPoptimal = cell(length(varFactors), simulationTime);
ratesMWBR = cell(length(varFactors), simulationTime);
ratesSCPF = cell(length(varFactors), simulationTime);
opBelongs = cell(length(varFactors), simulationTime);

pOutageSCPF = zeros(1, length(varFactors)); % it's an outage as long as there is one user not meeting minreq.
pOutageGPS = zeros(1, length(varFactors));
pOutageDP = zeros(1, length(varFactors));
pOutageDPWF = zeros(1, length(varFactors));
pOutageDPoptimal = zeros(1, length(varFactors));
pOutageMWBR = zeros(1, length(varFactors));

meanBtdGPS = zeros(1, length(varFactors));
meanBtdDP = zeros(1, length(varFactors));
meanBtdDPWF = zeros(1, length(varFactors));
meanBtdDPoptimal = zeros(1, length(varFactors));
meanBtdMWBR = zeros(1, length(varFactors));
meanBtdSCPF = zeros(1, length(varFactors));
meanUtilGPS = zeros(1, length(varFactors));
meanUtilDP = zeros(1, length(varFactors));
meanUtilDPWF = zeros(1, length(varFactors));
meanUtilDPoptimal = zeros(1, length(varFactors));
meanUtilMWBR = zeros(1, length(varFactors));
meanUtilSCPF = zeros(1, length(varFactors));

%% Run simulations
for i = 1:length(varFactors)
    varFactor = varFactors(i);
    rhoVec = relativeRhoVec * varFactor;
    bsAssociation = cell(1, simulationTime);
    capacities = cell(1, simulationTime);
    minRateReq = 0.4 * capacity / (varFactor * perBSLoad) * ones(1, nSlice); % min rate requirement
    shareDist = sharedimension(minRateReq, rhoVec, shareVec, outageTol, ...
        minSharePerBS, varFactor, 0);
    
    outageSCPF = zeros(1, simulationTime); 
    outageGPS = zeros(1, simulationTime);
    outageDP = zeros(1, simulationTime);
    outageDPWF = zeros(1, simulationTime);
    outageDPoptimal = zeros(1, simulationTime);
    outageMWBR = zeros(1, simulationTime);
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
        
        perUserMinRateReq = zeros(1, nUsers);
        for v = 1:nSlice
            perUserMinRateReq(opVec == v) = minRateReq(v);
        end
        
%         [r, f, b] = flexibleGPS(tmpNetSettings, tmpOpSettings, capacities{t}, ...
%             bsAssociation{t}, perUserMinRateReq);
        [r, f, b] = GPS(tmpNetSettings, tmpOpSettings, capacities{t}, bsAssociation{t});
        ratesGPS{i, t} = r;
        outageGPS(t) = any(r < unique(minRateReq));
        [r, f, b] = SCPF(tmpNetSettings, tmpOpSettings, capacities{t}, ...
            bsAssociation{t});
        ratesSCPF{i, t} = r;
        outageSCPF(t) = any(r < unique(minRateReq));
        [r, f, b] = DIFFPRICE(tmpNetSettings, tmpOpSettings, capacities{t}, ...
            bsAssociation{t}, minRateReq, 0);
        ratesDP{i, t} = r;
        outageDP(t) = any(r < unique(minRateReq));
        [r, f, b] = DIFFPRICE(tmpNetSettings, tmpOpSettings, capacities{t}, ...
            bsAssociation{t}, minRateReq, 1);
        ratesDPWF{i, t} = r;
        outageDPWF(t) = any(r < unique(minRateReq));
        if (sum(ratesDPWF{i, t} < 1e-4) > 0)
            ratesDPWF{i, t}(ratesDPWF{i, t} < 1e-4) = nan;
        end
        [r, f, b] = DPoptimal(tmpNetSettings, tmpOpSettings, capacities{t}, ...
            bsAssociation{t});
        ratesDPoptimal{i, t} = r;
        outageDPoptimal(t) = any(r < unique(minRateReq));
        [r, f, b] = MAXWEIGHT_BR(tmpNetSettings, tmpOpSettings, capacities{t}, ...
            bsAssociation{t});
        ratesMWBR{i, t} = r;
        outageMWBR(t) = any(r < unique(minRateReq));
        if (sum(ratesMWBR{i, t} < 1e-4) > 0)
            ratesMWBR{i, t}(ratesMWBR{i, t} < 1e-4) = nan;
        end
        fprintf('finish at time %d\n', t);
    end
    
    pOutageDP(i) = sum(outageDP) / simulationTime;
    pOutageDPWF(i) = sum(outageDPWF) / simulationTime;
    pOutageDPoptimal(i) = sum(outageDPoptimal) / simulationTime;
    pOutageGPS(i) = sum(outageGPS) / simulationTime;
    pOutageMWBR(i) = sum(outageMWBR) / simulationTime;
    pOutageSCPF(i) = sum(outageSCPF) / simulationTime;
    
    flatRateGPS = horzcat(ratesGPS{i, :});
    flatRateDPWF = horzcat(ratesDPWF{i, :});
    flatRateDP = horzcat(ratesDP{i, :});
    flatRateDPoptimal = horzcat(ratesDPoptimal{i, :});
    flatRateMWBR = horzcat(ratesMWBR{i, :});
    flatRateSCPF = horzcat(ratesSCPF{i, :});
    meanBtdGPS(i) = mean(1./flatRateGPS);
    meanBtdDP(i) = mean(1./flatRateDP);
    meanBtdDPWF(i) = mean(1./flatRateDPWF);
    meanBtdDPoptimal(i) = mean(1./flatRateDPoptimal);
    meanBtdMWBR(i) = nanmean(1./flatRateMWBR);
    meanBtdSCPF(i) = mean(1./flatRateSCPF);
  
    btdGainVecSCPF(i) = mean(1./flatRateGPS) / mean(1./flatRateSCPF); 
    btdGainVecDP(i) = mean(1./flatRateGPS) / mean(1./flatRateDP);
    btdGainVecDPWF(i) = mean(1./flatRateGPS) ...
        / nanmean(1./flatRateDPWF(flatRateDPWF > 1e-4));
    btdGainVecDPoptimal(i) = mean(1./flatRateGPS) / mean(1./flatRateDPoptimal);
    btdGainVecMWBR(i) = mean(1./flatRateGPS) ...
        / nanmean(1./flatRateMWBR(flatRateMWBR > 1e-4));
    utilGPS = zeros(1, simulationTime);
    utilSCPF = zeros(1, simulationTime);
    utilDP = zeros(1, simulationTime);
    utilDPWF = zeros(1, simulationTime);
    utilDPoptimal = zeros(1, simulationTime);
    utilMWBR = zeros(1, simulationTime);
    
    parfor t = 1:simulationTime % stats
        utilGPS(t) = ratetoutil(ratesGPS{i, t}, shareVec, opBelongs{i, t});
        utilSCPF(t) = ratetoutil(ratesSCPF{i, t}, shareVec, opBelongs{i, t});
        utilDP(t) = ratetoutil(ratesDP{i, t}, shareVec, opBelongs{i, t});
        if (sum(ratesDPWF{i, t} < 1e-4) > 0)
            utilDPWF(t) = nan;
        else
            utilDPWF(t) = ratetoutil(ratesDPWF{i, t}, shareVec, opBelongs{i, t});
        end
        utilDPoptimal(t) = ratetoutil(ratesDPoptimal{i, t}, shareVec, opBelongs{i, t});
        if (sum(ratesMWBR{i, t} < 1e-4) > 0)
            utilMWBR(t) = nan;
        else
            utilMWBR(t) = ratetoutil(ratesMWBR{i, t}, shareVec, opBelongs{i, t});
        end
    end
    meanUtilGPS(i) = nanmean(utilGPS);
    meanUtilDP(i) = nanmean(utilDP);
    meanUtilDPWF(i) = nanmean(utilDPWF);
    meanUtilDPoptimal(i) = nanmean(utilDPoptimal);
    meanUtilMWBR(i) = nanmean(utilMWBR);
    meanUtilSCPF(i) = nanmean(utilSCPF);
    utilityGainVecSCPF(i) = nanmean(utilGPS) / nanmean(utilSCPF); 
    utilityGainVecDP(i) = nanmean(utilGPS) / nanmean(utilDP);
    utilityGainVecDPWF(i) = nanmean(utilGPS) / nanmean(utilDPWF);
    utilityGainVecDPoptimal(i) = nanmean(utilGPS) / nanmean(utilDPoptimal);
    utilityGainVecMWBR(i) = nanmean(utilGPS) / nanmean(utilMWBR);
    
end

% Plot results
datestring = datestr(now, 30);

benchmarks = {'SCPF', 'DIFFPRICE-equal surplus', 'DIFFPRICE-waterfill', 'DIFFPRICE-optimal', ...
    'MAXWEIGHT-best response', 'GPS'};
bmWoGPS = {'SCPF', 'DIFFPRICE-equal surplus', 'DIFFPRICE-waterfill', 'DIFFPRICE-optimal', ...
    'MAXWEIGHT-best response'};

figure(7)
grid on
hold on
plot(varFactors, pOutageSCPF, 'b+-');
plot(varFactors, pOutageDP, 'ro-');
plot(varFactors, pOutageDPWF, 'bx:');
plot(varFactors, pOutageDPoptimal, 'ch-');
plot(varFactors, pOutageMWBR, 'kx-');
plot(varFactors, pOutageGPS, 'gd-');
title('P(outage) vs. mean load factor');
xlabel('mean load factor');
legend(benchmarks);
savefig(sprintf('figs/poutage-vs-var-%s.fig', datestring));

figure(1);
grid on
hold on
plot(varFactors, btdGainVecSCPF, 'b+-');
plot(varFactors, btdGainVecDP, 'ro-');
plot(varFactors, btdGainVecDPWF, 'bx:');
plot(varFactors, btdGainVecDPoptimal, 'ch-');
plot(varFactors, btdGainVecMWBR, 'kx-');
title('BTD gain over GPS vs. mean load factor');
xlabel('mean load factor');
legend(bmWoGPS);
savefig(sprintf('figs/btd-gain-vs-var-%s.fig', datestring));

figure(2);
grid on
hold on
plot(varFactors, meanUtilSCPF - meanUtilGPS, 'b+-');
plot(varFactors, meanUtilDP - meanUtilGPS, 'ro-');
plot(varFactors, meanUtilDPWF - meanUtilGPS, 'bx:');
plot(varFactors, meanUtilDPoptimal - meanUtilGPS, 'ch-');
plot(varFactors, meanUtilMWBR - meanUtilGPS, 'kx-');
title('Utility gain over GPS vs. mean load factor');
xlabel('mean load factor');
legend(bmWoGPS);
savefig(sprintf('figs/util-gain-vs-var-%s.fig', datestring));

% Get some idea on slice idx
idx = 2;
btdGainVecSCPF1 = zeros(1, length(varFactors)); % BTD gain over (flexible) GPS.
btdGainVecDP1 = zeros(1, length(varFactors));
btdGainVecDPWF1 = zeros(1, length(varFactors));
btdGainVecDPoptimal1 = zeros(1, length(varFactors));
btdGainVecMWBR1 = zeros(1, length(varFactors));
meanBtdGPS1 = zeros(1, length(varFactors));
meanBtdDP1 = zeros(1, length(varFactors));
meanBtdDPWF1 = zeros(1, length(varFactors));
meanBtdDPoptimal1 = zeros(1, length(varFactors));
meanBtdMWBR1 = zeros(1, length(varFactors));
meanBtdSCPF1 = zeros(1, length(varFactors));

for i = 1:length(varFactors)
    sliceIdx = (horzcat(opBelongs{i, :}) == idx);
    flatRateGPS1 = horzcat(ratesGPS{i, :});
    flatRateGPS1 = flatRateGPS1(sliceIdx);
    flatRateDP1 = horzcat(ratesDP{i, :});
    flatRateDP1 = flatRateDP1(sliceIdx);
    flatRateDPWF1 = horzcat(ratesDPWF{i, :});
    flatRateDPWF1 = flatRateDPWF1(sliceIdx);
    flatRateDPoptimal1 = horzcat(ratesDPoptimal{i, :});
    flatRateDPoptimal1 = flatRateDPoptimal1(sliceIdx);
    flatRateMWBR1 = horzcat(ratesMWBR{i, :});
    flatRateMWBR1 = flatRateMWBR1(sliceIdx);
    flatRateMWBR1 = flatRateMWBR1(flatRateMWBR1 > 1e-4);
    flatRateSCPF1 = horzcat(ratesSCPF{i, :});
    flatRateSCPF1 = flatRateSCPF1(sliceIdx);
    meanBtdGPS1(i) = mean(1./flatRateGPS1);
    meanBtdDP1(i) = mean(1./flatRateDP1);
    meanBtdDPWF1(i) = mean(1./flatRateDPWF1);
    meanBtdDPoptimal1(i) = mean(1./flatRateDPoptimal1);
    meanBtdMWBR1(i) = nanmean(1./flatRateMWBR1);
    meanBtdSCPF1(i) = mean(1./flatRateSCPF1);
    
    btdGainVecSCPF1(i) = mean(1./flatRateGPS1) / mean(1./flatRateSCPF1); 
    btdGainVecDP1(i) = mean(1./flatRateGPS1) / mean(1./flatRateDP1);
    btdGainVecDPWF1(i) = mean(1./flatRateGPS1) / nanmean(1./flatRateDPWF1);
    btdGainVecDPoptimal1(i) = mean(1./flatRateGPS1) / mean(1./flatRateDPoptimal1);
    btdGainVecMWBR1(i) = mean(1./flatRateGPS1) ...
        / nanmean(1./flatRateMWBR1);
end

figure(3)
grid on
hold on
plot(varFactors, meanBtdSCPF1, 'b+-');
plot(varFactors, meanBtdDP1, 'ro-');
plot(varFactors, meanBtdDPWF1, 'bx:');
plot(varFactors, meanBtdDPoptimal1, 'ch-');
plot(varFactors, meanBtdMWBR1, 'kx-');
plot(varFactors, meanBtdGPS1, 'gd-');
title('Average btd vs. mean load factor of slice 2');
xlabel('mean load factor');
legend(benchmarks);
savefig(sprintf('figs/btd-vs-var-slice2-%s.fig', datestring));

figure(4)
grid on
hold on
plot(varFactors, btdGainVecSCPF1, 'b+-');
plot(varFactors, btdGainVecDP1, 'ro-');
plot(varFactors, btdGainVecDPWF1, 'bx:');
plot(varFactors, btdGainVecDPoptimal1, 'ch-');
plot(varFactors, btdGainVecMWBR1, 'kx-');
title('BTD gain over GPS vs. mean load factor of slice 2');
xlabel('mean load factor');
legend(bmWoGPS);
savefig(sprintf('figs/btd-gain-vs-var-slice2-%s.fig', datestring));

idx = 1;
btdGainVecSCPF1 = zeros(1, length(varFactors)); % BTD gain over (flexible) GPS.
btdGainVecDP1 = zeros(1, length(varFactors));
btdGainVecDPWF1 = zeros(1, length(varFactors));
btdGainVecDPoptimal1 = zeros(1, length(varFactors));
btdGainVecMWBR1 = zeros(1, length(varFactors));
meanBtdGPS1 = zeros(1, length(varFactors));
meanBtdDP1 = zeros(1, length(varFactors));
meanBtdDPWF1 = zeros(1, length(varFactors));
meanBtdDPoptimal1 = zeros(1, length(varFactors));
meanBtdMWBR1 = zeros(1, length(varFactors));
meanBtdSCPF1 = zeros(1, length(varFactors));

for i = 1:length(varFactors)
    sliceIdx = (horzcat(opBelongs{i, :}) == idx);
    flatRateGPS1 = horzcat(ratesGPS{i, :});
    flatRateGPS1 = flatRateGPS1(sliceIdx);
    flatRateDP1 = horzcat(ratesDP{i, :});
    flatRateDP1 = flatRateDP1(sliceIdx);
    flatRateDPWF1 = horzcat(ratesDPWF{i, :});
    flatRateDPWF1 = flatRateDPWF1(sliceIdx);
    flatRateDPoptimal1 = horzcat(ratesDPoptimal{i, :});
    flatRateDPoptimal1 = flatRateDPoptimal1(sliceIdx);
    flatRateMWBR1 = horzcat(ratesMWBR{i, :});
    flatRateMWBR1 = flatRateMWBR1(sliceIdx);
    flatRateMWBR1 = flatRateMWBR1(flatRateMWBR1 > 1e-4);
    flatRateSCPF1 = horzcat(ratesSCPF{i, :});
    flatRateSCPF1 = flatRateSCPF1(sliceIdx);
    meanBtdGPS1(i) = mean(1./flatRateGPS1);
    meanBtdDP1(i) = mean(1./flatRateDP1);
    meanBtdDPWF1(i) = mean(1./flatRateDPWF1);
    meanBtdDPoptimal1(i) = mean(1./flatRateDPoptimal1);
    meanBtdMWBR1(i) = nanmean(1./flatRateMWBR1);
    meanBtdSCPF1(i) = mean(1./flatRateSCPF1);
    
    btdGainVecSCPF1(i) = mean(1./flatRateGPS1) / mean(1./flatRateSCPF1); 
    btdGainVecDP1(i) = mean(1./flatRateGPS1) / mean(1./flatRateDP1);
    btdGainVecDPWF1(i) = mean(1./flatRateGPS1) / nanmean(1./flatRateDPWF1);
    btdGainVecDPoptimal1(i) = mean(1./flatRateGPS1) / mean(1./flatRateDPoptimal1);
    btdGainVecMWBR1(i) = mean(1./flatRateGPS1) ...
        / nanmean(1./flatRateMWBR1);
end

figure(5)
grid on
hold on
plot(varFactors, meanBtdSCPF1, 'b+-');
plot(varFactors, meanBtdDP1, 'ro-');
plot(varFactors, meanBtdDPWF1, 'bx:');
plot(varFactors, meanBtdDPoptimal1, 'ch-');
plot(varFactors, meanBtdMWBR1, 'kx-');
plot(varFactors, meanBtdGPS1, 'gd-');
title('Average btd vs. mean load factor of slice 1');
xlabel('mean load factor');
legend(benchmarks);
savefig(sprintf('figs/btd-vs-var-slice1-%s.fig', datestring));

figure(6)
grid on
hold on
plot(varFactors, btdGainVecSCPF1, 'b+-');
plot(varFactors, btdGainVecDP1, 'ro-');
plot(varFactors, btdGainVecDPWF1, 'bx:');
plot(varFactors, btdGainVecDPoptimal1, 'ch-');
plot(varFactors, btdGainVecMWBR1, 'kx-');
title('BTD gain over GPS vs. mean load factor of slice 1');
xlabel('mean load factor');
legend(bmWoGPS);
savefig(sprintf('figs/btd-gain-vs-var-slice1-%s.fig', datestring));
