clc, close all, clear all
nSlice = 4;
%parpool('local', 40);
warning('off','all');

simulationTime = 4;
% Setup:
% Two inelastic slices, with uniform minimal rate requirements, and no 
% elasticity in the utility function.
% Another two are elastic users, with 0 minimal rate requirements,
% elasticity is specified as 1-fairness, i.e., log utility. 
% We simulate the scenario where Slices 1 and 2 are congested at BSs 1 and
% 2, where elastic users are roughly uniform.
perBSLoad = 1;
% shareVec = [1 1 1 1];
sliceCats = [0 0 1 1];
afterSliceCats = [1 1 1 1];
relativeRhoVec = perBSLoad * [[2 2 6 6];
                              [10 10 6 6];
                              [10 10 6 6];
                              [2 2 6 6]]';
nBaseStations = size(relativeRhoVec, 2);
capacity = 1;
minSharePerBS = 0.001;
outageTol = 0.01;
netSettings = [];
netSettings.bsNS = nBaseStations;
opSettings = [];

meanFactorVec = 1:0.5:4;

btdGainVecSCPF = zeros(1, length(meanFactorVec)); % BTD gain over (flexible) GPS.
btdGainVecDP = zeros(1, length(meanFactorVec));
btdGainVecDPoptimal = zeros(1, length(meanFactorVec));

utilityGainVecSCPF = zeros(1, length(meanFactorVec)); % overall utility gain over (flexible) GPS.
utilityGainVecDP = zeros(1, length(meanFactorVec));
utilityGainVecDPoptimal = zeros(1, length(meanFactorVec));

ratesGPS = cell(length(meanFactorVec), simulationTime); % Save ordinary data for regression.
ratesDP = cell(length(meanFactorVec), simulationTime);
ratesDPoptimal = cell(length(meanFactorVec), simulationTime);
ratesSCPF = cell(length(meanFactorVec), simulationTime);

weightsPerUsersMat = cell(length(meanFactorVec), simulationTime);

opBelongs = cell(length(meanFactorVec), simulationTime);

pOutageSCPF = zeros(1, length(meanFactorVec)); % it's an outage as long as there is one user not meeting minreq.
pOutageGPS = zeros(1, length(meanFactorVec));
pOutageDP = zeros(1, length(meanFactorVec));
pOutageDPoptimal = zeros(1, length(meanFactorVec));

meanBtdGPS = zeros(1, length(meanFactorVec));
meanBtdDP = zeros(1, length(meanFactorVec));
meanBtdDPoptimal = zeros(1, length(meanFactorVec));
meanBtdSCPF = zeros(1, length(meanFactorVec));

meanUtilGPS = zeros(1, length(meanFactorVec));
meanUtilDP = zeros(1, length(meanFactorVec));
meanUtilDPoptimal = zeros(1, length(meanFactorVec));
meanUtilSCPF = zeros(1, length(meanFactorVec));

for i = 1:length(meanFactorVec)
    varFactor = 1;
    rhoVec = relativeRhoVec * varFactor;
    bsAssociation = cell(1, simulationTime);
    capacities = cell(1, simulationTime);
    minRateReq = 0.025 * capacity / (varFactor * perBSLoad) * ones(1, nSlice);
    minRateReq(3:4) = meanFactorVec(i) * minRateReq(3:4);
    afterMinRate = minRateReq;
    afterMinRate(sliceCats > 0) = 0;
    
    [shareDist, gpsShareDist, shareVec] = sharedimension(minRateReq, rhoVec, outageTol, ...
            minSharePerBS, 1, 0, sliceCats, ones(1, 4), ones(4, 4));
    
    outageSCPF = zeros(1, simulationTime); 
    outageGPS = zeros(1, simulationTime);
    outageDP = zeros(1, simulationTime);
    outageDPoptimal = zeros(1, simulationTime);
    
    totalNumUsers = 0;
    
    parfor t = 1:simulationTime
        loadDist = poissrnd(rhoVec);
        nUsers = sum(sum(loadDist));
        totalNumUsers = totalNumUsers + nUsers;
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
        tmpOpSettings.s_o = shareVec;
        tmpOpSettings.ops_belongs = opBelongs{i, t};
        tmpOpSettings.shareDist = shareDist;
        weightPerUser = zeros(1, tmpNetSettings.users);
        for v = 1:nSlice
            weightPerUser(tmpOpSettings.ops_belongs == v) = shareVec(v) ...
                ./ sum(tmpOpSettings.ops_belongs == v);
        end
        tmpOpSettings.w_i = weightPerUser;
        weightsPerUsersMat{i, t} = weightPerUser;
        
        perUserMinRateReq = zeros(1, nUsers);
        for v = 1:nSlice
            if (sliceCats(v) == 0)
                perUserMinRateReq(tmpOpSettings.ops_belongs == v) = minRateReq(v);
            end
        end
        
        
        
        [r, f, b] = SCPF(tmpNetSettings, tmpOpSettings, capacities{t}, ...
            bsAssociation{t});
        ratesSCPF{i, t} = r;
        outageSCPF(t) = sum(r < perUserMinRateReq);
        
        [r, f, b] = DIFFPRICE(tmpNetSettings, tmpOpSettings, capacities{t}, ...
            bsAssociation{t}, afterMinRate, 0, afterSliceCats, weightPerUser);
        ratesDP{i, t} = r;
        outageDP(t) = sum(r < perUserMinRateReq);
        
        [r, f, b] = DPoptimal(tmpNetSettings, tmpOpSettings, capacities{t}, ...
            bsAssociation{t}, perUserMinRateReq, afterSliceCats);
        ratesDPoptimal{i, t} = r;
        outageDPoptimal(t) = sum(r < perUserMinRateReq);
        
        % GPS, needs to first adjust the share dimensioning.
        tmpOpSettings.shareDist = gpsShareDist;
        [r, f, b] = flexibleGPS(tmpNetSettings, tmpOpSettings, capacities{t}, ...
            bsAssociation{t}, perUserMinRateReq); % dummy minreq.
        ratesGPS{i, t} = r;
        outageGPS(t) = sum(r < perUserMinRateReq);
        
        fprintf('finish at time %d\n', t);
    end
    
    pOutageDP(i) = sum(outageDP) / totalNumUsers;
    pOutageDPoptimal(i) = sum(outageDPoptimal) / totalNumUsers;
    pOutageGPS(i) = sum(outageGPS) / totalNumUsers;
    pOutageSCPF(i) = sum(outageSCPF) / totalNumUsers;
    
    flatRateGPS = horzcat(ratesGPS{i, :});
    flatRateDP = horzcat(ratesDP{i, :});
    flatRateDPoptimal = horzcat(ratesDPoptimal{i, :});
    flatRateSCPF = horzcat(ratesSCPF{i, :});
    
    meanBtdGPS(i) = mean(1./flatRateGPS);
    meanBtdDP(i) = mean(1./flatRateDP);
    meanBtdDPoptimal(i) = mean(1./flatRateDPoptimal);
    meanBtdSCPF(i) = mean(1./flatRateSCPF);
  
    btdGainVecSCPF(i) = mean(1./flatRateGPS) / mean(1./flatRateSCPF); 
    btdGainVecDP(i) = mean(1./flatRateGPS) / mean(1./flatRateDP);
    btdGainVecDPoptimal(i) = mean(1./flatRateGPS) / mean(1./flatRateDPoptimal);
    
    utilGPS = zeros(1, simulationTime);
    utilSCPF = zeros(1, simulationTime);
    utilDP = zeros(1, simulationTime);
    utilDPoptimal = zeros(1, simulationTime);
    
    
    for t = 1:simulationTime % stats
        % For each time instant, only account for the set of users receiving at
        % least min rate req under all benchmarks.
        nUsers = length(capacities{t});
        opVec = opBelongs{i, t};
        perUserMinRateReq = zeros(1, nUsers);
        for v = 1:nSlice
            perUserMinRateReq(opVec == v) = afterMinRate(v);
        end
        
        goodUsers = (ratesGPS{i, t} > perUserMinRateReq & ratesSCPF{i, t} ...
            > perUserMinRateReq & ratesDP{i, t} > perUserMinRateReq ...
            & ratesDPoptimal{i, t} > perUserMinRateReq);
        
        tmpRatesGPS = nan(size(ratesGPS{i, t}));
        tmpRatesGPS(goodUsers) = ratesGPS{i, t}(goodUsers);
        utilGPS(t) = ratetoutil_old(tmpRatesGPS, shareVec, ...
            opBelongs{i, t}, afterSliceCats, perUserMinRateReq, weightsPerUsersMat{i, t});
        
        tmpRatesSCPF = nan(size(ratesSCPF{i, t}));
        tmpRatesSCPF(goodUsers) = ratesSCPF{i, t}(goodUsers);
        utilSCPF(t) = ratetoutil_old(tmpRatesSCPF, shareVec, ...
            opBelongs{i, t}, afterSliceCats, perUserMinRateReq, weightsPerUsersMat{i, t});
        
        tmpRatesDP = nan(size(ratesDP{i, t}));
        tmpRatesDP(goodUsers) = ratesDP{i, t}(goodUsers);
        utilDP(t) = ratetoutil_old(tmpRatesDP, shareVec, ...
            opBelongs{i, t}, afterSliceCats, perUserMinRateReq, weightsPerUsersMat{i, t});
        
        tmpRatesDPoptimal = nan(size(ratesDPoptimal{i, t}));
        tmpRatesDPoptimal(goodUsers) = ratesDPoptimal{i, t}(goodUsers);
        utilDPoptimal(t) = ratetoutil_old(tmpRatesDPoptimal, shareVec, ...
            opBelongs{i, t}, afterSliceCats, perUserMinRateReq, weightsPerUsersMat{i, t});
        
%         if (outageGPS(t))
%             utilGPS(t) = nan;
%         else
%             utilGPS(t) = ratetoutil(ratesGPS{i, t}, shareVec, ...
%                 opBelongs{i, t}, sliceCats, perUserMinRateReq);
%         end
%         if (outageSCPF(t))
%             utilSCPF(t) = nan;
%         else
%             utilSCPF(t) = ratetoutil(ratesSCPF{i, t}, shareVec, ...
%                 opBelongs{i, t}, sliceCats, perUserMinRateReq);
%         end
%         if (outageDP(t))
%             utilDP(t) = nan;
%         else
%             utilDP(t) = ratetoutil(ratesDP{i, t}, shareVec, opBelongs{i, t}, ...
%                 sliceCats, perUserMinRateReq);
%         end
%         if (outageDPoptimal(t))
%             utilDPoptimal(t) = nan;
%         else
%             utilDPoptimal(t) = ratetoutil(ratesDPoptimal{i, t}, shareVec, ...
%                 opBelongs{i, t}, sliceCats, perUserMinRateReq);
%         end
    end
    
    meanUtilGPS(i) = nanmean(utilGPS);
    meanUtilDP(i) = nanmean(utilDP);
    meanUtilDPoptimal(i) = nanmean(utilDPoptimal);
    meanUtilSCPF(i) = nanmean(utilSCPF);
    
end

%% Plot something
datestring = datestr(now, 30);

benchmarks = {'Share-based approach', 'GREET', 'Social optimal', ...
    'Reservation-based approach'};
bmWoGPS = {'Share-based approach', 'GREET', 'Social optimal'};

figure(8)
grid on
hold on
h1 = plot(meanUtilSCPF, 1 - pOutageSCPF, 'b+-');
h2 = plot(meanUtilDP, 1 - pOutageDP, 'ro-');
h3 = plot(meanUtilDPoptimal, 1 - pOutageDPoptimal, 'ch-');
h4 = plot(meanUtilGPS, 1 - pOutageGPS, 'gd-');
set([h1 h2 h3 h4], 'LineWidth', 2, 'MarkerSize', 10);
title('P(outage) vs. utility trade off');
xlabel('utility');
ylabel('1 - P(outage)')
legend(benchmarks);
savefig(sprintf('figs/util-outage-tradeoff-%s.fig', datestring));

% figure(7)
% grid on
% hold on
% plot(varFactors, pOutageSCPF, 'b+-');
% plot(varFactors, pOutageDP, 'ro-');
% plot(varFactors, pOutageDPoptimal, 'ch-');
% plot(varFactors, pOutageGPS, 'gd-');
% title('P(outage) vs. mean load factor');
% xlabel('mean load factor');
% legend(benchmarks);
% savefig(sprintf('figs/poutage-vs-var-%s.fig', datestring));
% 
% figure(1);
% grid on
% hold on
% plot(varFactors, btdGainVecSCPF, 'b+-');
% plot(varFactors, btdGainVecDP, 'ro-');
% plot(varFactors, btdGainVecDPoptimal, 'ch-');
% title('BTD gain over GPS vs. mean load factor');
% xlabel('mean load factor');
% legend(bmWoGPS);
% savefig(sprintf('figs/btd-gain-vs-var-%s.fig', datestring));
% 
% figure(2);
% grid on
% hold on
% plot(varFactors, meanUtilSCPF - meanUtilGPS, 'b+-');
% plot(varFactors, meanUtilDP - meanUtilGPS, 'ro-');
% plot(varFactors, meanUtilDPoptimal - meanUtilGPS, 'ch-');
% title('Utility gain over GPS vs. mean load factor');
% xlabel('mean load factor');
% legend(bmWoGPS);
% savefig(sprintf('figs/util-gain-vs-var-%s.fig', datestring));
% 
% % Get some idea on slice idx
% idx = 2;
% 
% btdGainVecSCPF1 = zeros(1, length(varFactors)); % BTD gain over (flexible) GPS.
% btdGainVecDP1 = zeros(1, length(varFactors));
% btdGainVecDPoptimal1 = zeros(1, length(varFactors));
% 
% meanBtdGPS1 = zeros(1, length(varFactors));
% meanBtdDP1 = zeros(1, length(varFactors));
% meanBtdDPoptimal1 = zeros(1, length(varFactors));
% meanBtdSCPF1 = zeros(1, length(varFactors));
% 
% for i = 1:length(varFactors)
%     sliceIdx = (horzcat(opBelongs{i, :}) == idx);
%     flatRateGPS1 = horzcat(ratesGPS{i, :});
%     flatRateGPS1 = flatRateGPS1(sliceIdx);
%     flatRateDP1 = horzcat(ratesDP{i, :});
%     flatRateDP1 = flatRateDP1(sliceIdx);
%     flatRateDPoptimal1 = horzcat(ratesDPoptimal{i, :});
%     flatRateDPoptimal1 = flatRateDPoptimal1(sliceIdx);
%     flatRateSCPF1 = horzcat(ratesSCPF{i, :});
%     flatRateSCPF1 = flatRateSCPF1(sliceIdx);
%     
%     meanBtdGPS1(i) = mean(1./flatRateGPS1);
%     meanBtdDP1(i) = mean(1./flatRateDP1); 
%     meanBtdDPoptimal1(i) = mean(1./flatRateDPoptimal1);
%     meanBtdSCPF1(i) = mean(1./flatRateSCPF1);
%     
%     btdGainVecSCPF1(i) = mean(1./flatRateGPS1) / mean(1./flatRateSCPF1); 
%     btdGainVecDP1(i) = mean(1./flatRateGPS1) / mean(1./flatRateDP1);
%     btdGainVecDPoptimal1(i) = mean(1./flatRateGPS1) / mean(1./flatRateDPoptimal1);
% 
% end
% 
% figure(3)
% grid on
% hold on
% plot(varFactors, meanBtdSCPF1, 'b+-');
% plot(varFactors, meanBtdDP1, 'ro-');
% plot(varFactors, meanBtdDPoptimal1, 'ch-');
% plot(varFactors, meanBtdGPS1, 'gd-');
% title('Average btd vs. mean load factor of slice 2');
% xlabel('mean load factor');
% legend(benchmarks);
% savefig(sprintf('figs/btd-vs-var-slice2-%s.fig', datestring));
% 
% figure(4)
% grid on
% hold on
% plot(varFactors, btdGainVecSCPF1, 'b+-');
% plot(varFactors, btdGainVecDP1, 'ro-');
% plot(varFactors, btdGainVecDPoptimal1, 'ch-');
% title('BTD gain over GPS vs. mean load factor of slice 2');
% xlabel('mean load factor');
% legend(bmWoGPS);
% savefig(sprintf('figs/btd-gain-vs-var-slice2-%s.fig', datestring));
% 
% idx = 3;
% btdGainVecSCPF1 = zeros(1, length(varFactors)); % BTD gain over (flexible) GPS.
% btdGainVecDP1 = zeros(1, length(varFactors));
% btdGainVecDPoptimal1 = zeros(1, length(varFactors));
% 
% meanBtdGPS1 = zeros(1, length(varFactors));
% meanBtdDP1 = zeros(1, length(varFactors));
% meanBtdDPoptimal1 = zeros(1, length(varFactors));
% meanBtdSCPF1 = zeros(1, length(varFactors));
% 
% for i = 1:length(varFactors)
%     sliceIdx = (horzcat(opBelongs{i, :}) == idx);
%     flatRateGPS1 = horzcat(ratesGPS{i, :});
%     flatRateGPS1 = flatRateGPS1(sliceIdx);
%     flatRateDP1 = horzcat(ratesDP{i, :});
%     flatRateDP1 = flatRateDP1(sliceIdx);
%     flatRateDPoptimal1 = horzcat(ratesDPoptimal{i, :});
%     flatRateDPoptimal1 = flatRateDPoptimal1(sliceIdx);
%     flatRateSCPF1 = horzcat(ratesSCPF{i, :});
%     flatRateSCPF1 = flatRateSCPF1(sliceIdx);
%     
%     meanBtdGPS1(i) = mean(1./flatRateGPS1);
%     meanBtdDP1(i) = mean(1./flatRateDP1);
%     meanBtdDPoptimal1(i) = mean(1./flatRateDPoptimal1);
%     meanBtdSCPF1(i) = mean(1./flatRateSCPF1);
%     
%     btdGainVecSCPF1(i) = mean(1./flatRateGPS1) / mean(1./flatRateSCPF1); 
%     btdGainVecDP1(i) = mean(1./flatRateGPS1) / mean(1./flatRateDP1);
%     btdGainVecDPoptimal1(i) = mean(1./flatRateGPS1) / mean(1./flatRateDPoptimal1);
% 
% end
% 
% figure(5)
% grid on
% hold on
% plot(varFactors, meanBtdSCPF1, 'b+-');
% plot(varFactors, meanBtdDP1, 'ro-');
% plot(varFactors, meanBtdDPoptimal1, 'ch-');
% plot(varFactors, meanBtdGPS1, 'gd-');
% title('Average btd vs. mean load factor of slice 3');
% xlabel('mean load factor');
% legend(benchmarks);
% savefig(sprintf('figs/btd-vs-var-slice1-%s.fig', datestring));
% 
% figure(6)
% grid on
% hold on
% plot(varFactors, btdGainVecSCPF1, 'b+-');
% plot(varFactors, btdGainVecDP1, 'ro-');
% plot(varFactors, btdGainVecDPoptimal1, 'ch-');
% title('BTD gain over GPS vs. mean load factor of slice 3');
% xlabel('mean load factor');
% legend(bmWoGPS);
% savefig(sprintf('figs/btd-gain-vs-var-slice1-%s.fig', datestring));
