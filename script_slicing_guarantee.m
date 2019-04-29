clc, close all, clear all
nSlice = 4;
parpool('local', 40);
warning('off','all');
nSlice = 3;

simulationTime = 4;
% Setup:
% Two inelastic slices, with uniform minimal rate requirements, and no 
% elasticity in the utility function.
% Another two are elastic users, with 0 minimal rate requirements,
% elasticity is specified as 1-fairness, i.e., log utility. 
% We simulate the scenario where Slices 1 and 2 are congested at BSs 1 and
% 2, where elastic users are roughly uniform.
perBSLoad = 6;
shareVec = [1 1 1 1];
sliceCats = [0 0 1 1];
relativeRhoVec = perBSLoad * [[2 2 6 6];
                              [10 10 6 6];
                              [10 10 6 6];
                              [2 2 6 6]]';
nBaseStations = size(relativeRhoVec, 2);
capacity = 1;
minSharePerBS = 0.01;
outageTol = 0.05;
netSettings = [];
netSettings.bsNS = nBaseStations;
opSettings = [];
opSettings.s_o = shareVec;

varFactors = 1:0.5:3;

btdGainVecSCPF = zeros(1, length(varFactors)); % BTD gain over (flexible) GPS.
btdGainVecDP = zeros(1, length(varFactors));
btdGainVecDPoptimal = zeros(1, length(varFactors));

utilityGainVecSCPF = zeros(1, length(varFactors)); % overall utility gain over (flexible) GPS.
utilityGainVecDP = zeros(1, length(varFactors));
utilityGainVecDPoptimal = zeros(1, length(varFactors));

ratesGPS = cell(length(varFactors), simulationTime); % Save ordinary data for regression.
ratesDP = cell(length(varFactors), simulationTime);
ratesDPoptimal = cell(length(varFactors), simulationTime);
ratesSCPF = cell(length(varFactors), simulationTime);

opBelongs = cell(length(varFactors), simulationTime);

pOutageSCPF = zeros(1, length(varFactors)); % it's an outage as long as there is one user not meeting minreq.
pOutageGPS = zeros(1, length(varFactors));
pOutageDP = zeros(1, length(varFactors));
pOutageDPoptimal = zeros(1, length(varFactors));

meanBtdGPS = zeros(1, length(varFactors));
meanBtdDP = zeros(1, length(varFactors));
meanBtdDPoptimal = zeros(1, length(varFactors));
meanBtdSCPF = zeros(1, length(varFactors));

meanUtilGPS = zeros(1, length(varFactors));
meanUtilDP = zeros(1, length(varFactors));
meanUtilDPoptimal = zeros(1, length(varFactors));
meanUtilSCPF = zeros(1, length(varFactors));

for i = 1:length(varFactors)
    rhoVec = relativeRhoVec * varFactor;
    bsAssociation = cell(1, simulationTime);
    capacities = cell(1, simulationTime);
    minRateReq = 0.4 * capacity / (varFactor * perBSLoad) * ones(1, nSlice);
    minRateReq(3:4) = 0;
    shareDist = sharedimension(minRateReq, rhoVec, shareVec, outageTol, ...
        minSharePerBS, varFactor, 0, 0);
    
    outageSCPF = zeros(1, simulationTime); 
    outageGPS = zeros(1, simulationTime);
    outageDP = zeros(1, simulationTime);
    outageDPoptimal = zeros(1, simulationTime);
    
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
        
        [r, f, b] = SCPF(tmpNetSettings, tmpOpSettings, capacities{t}, ...
            bsAssociation{t});
        ratesSCPF{i, t} = r;
        outageSCPF(t) = any(r < perUserMinRateReq);
        
        [r, f, b] = DIFFPRICE(tmpNetSettings, tmpOpSettings, capacities{t}, ...
            bsAssociation{t}, minRateReq, 0);
        ratesDP{i, t} = r;
        outageDP(t) = any(r < perUserMinRateReq);
        
        [r, f, b] = DPoptimal(tmpNetSettings, tmpOpSettings, capacities{t}, ...
            bsAssociation{t}, perUserMinRateReq, sliceCats);
        ratesDPoptimal{i, t} = r;
        outageDPoptimal(t) = any(r < perUserMinRateReq);
        
        % GPS, needs to first adjust the share dimensioning.
        tmpOpSettings.shareDist = sharedimension(minRateReq, rhoVec, ...
            shareVec, outageTol, minSharePerBS, varFactor, 0, 1);
        [r, f, b] = GPS(tmpNetSettings, tmpOpSettings, capacities{t}, bsAssociation{t});
        ratesGPS{i, t} = r;
        outageGPS(t) = any(r < perUserMinRateReq);
        
        fprintf('finish at time %d\n', t);
    end
    
    pOutageDP(i) = sum(outageDP) / simulationTime;
    pOutageDPoptimal(i) = sum(outageDPoptimal) / simulationTime;
    pOutageGPS(i) = sum(outageGPS) / simulationTime;
    pOutageSCPF(i) = sum(outageSCPF) / simulationTime;
    
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
    
    parfor t = 1:simulationTime % stats
        nUsers = length(capacities{t});
        opVec = opBelongs{i, t};
        perUserMinRateReq = zeros(1, nUsers);
        for v = 1:nSlice
            perUserMinRateReq(opVec == v) = minRateReq(v);
        end
        utilGPS(t) = ratetoutil(ratesGPS{i, t}, shareVec, ...
            opBelongs{i, t}, sliceCats, perUserMinRateReq);
        utilSCPF(t) = ratetoutil(ratesSCPF{i, t}, shareVec, ...
            opBelongs{i, t}, sliceCats, perUserMinRateReq);
        utilDP(t) = ratetoutil(ratesDP{i, t}, shareVec, opBelongs{i, t}, ...
            sliceCats, perUserMinRateReq);
        utilDPoptimal(t) = ratetoutil(ratesDPoptimal{i, t}, shareVec, ...
            opBelongs{i, t}, sliceCats, perUserMinRateReq);
    end
    
    meanUtilGPS(i) = nanmean(utilGPS);
    meanUtilDP(i) = nanmean(utilDP);
    meanUtilDPoptimal(i) = nanmean(utilDPoptimal);
    meanUtilSCPF(i) = nanmean(utilSCPF);
    
    utilityGainVecSCPF(i) = nanmean(utilGPS) / nanmean(utilSCPF); 
    utilityGainVecDP(i) = nanmean(utilGPS) / nanmean(utilDP);
    utilityGainVecDPoptimal(i) = nanmean(utilGPS) / nanmean(utilDPoptimal);
end