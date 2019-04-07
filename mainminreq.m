% mainminreq The script used to simulate the settings where minimal rate
% requirements are specified for users.
% Author: numerology
% 4/4/2019
clc, close all, clear all, gcp;
%% Settings
nSlice = 3;

simulationTime = 10;
perBSLoad = 6;
shareVec = [14/9 13/18 13/18];
relativeRhoVec = [perBSLoad * [2/3 1/6 1/6];
                  perBSLoad * [2/3 1/6 1/6];
                  perBSLoad * [2/9 7/18 7/18]]'; % mean load distribution, V x B
% Share dimensioning

shareDist = [2/3 1/6 1/6;
    2/3 1/6 1/6;
    2/9 7/18 7/18]';

nBaseStations = size(relativeRhoVec, 2);
capacity = 1;
minRateReq = 0.5 * capacity / perBSLoad * ones(1, nSlice); % min rate requirement
outageTol = 0.2;
netSettings = [];
netSettings.bsNS = nBaseStations;
opSettings = [];
opSettings.s_o = shareVec;
opSettings.shareDist = shareDist;

% share dimensioning
shareDist = sharedimension(minRateReq, relativeRhoVec, shareVec, outageTol);

%% Generate the user association
bsAssociation = cell(1, simulationTime);
opBelongs = cell(1, simulationTime);
capacities = cell(1, simulationTime);
ppm = ParforProgMon('Generating network profile: ', simulationTime);
parfor t = 1:simulationTime
    loadDist = 1 * poissrnd(relativeRhoVec);
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
ratesDP = cell(1, simulationTime);
ppm = ParforProgMon('Simulating resource sharing : ', simulationTime);
for t = 1:simulationTime
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
    [r, f, b] = DIFFPRICE(tmpNetSettings, tmpOpSettings, capacities{t}, ...
        bsAssociation{t}, minRateReq);
    ratesDP{t} = r;
    ppm.increment();
end

