function [NetSettings, OpSettings, c_uD, bsD, trace, bsPositions, ...
    userDemands]=Network_configuration(simulationTime, warmup, bsN, sectors,...
    interdistance, shareVec, phiLevels, sat, operators, alphas, ~)
% Generate the simulation scenarios
% simulationTime: duration for the simulation
% warmup: if begin with a warm start for link estimation, this is the
% starting time
% bsN: number of base stations
% sectors: each base station is split into how many sectors.
% interdistance: distance between base stations
% shareVec: vector of share allocation among slices 
% phiLevels: parameter used in intraslice user priority allocation
% sat: (avg) number of users per sector
% operators: number of slices
% alphas: alpha in the alpha-fairness notation, used in computing utility
% functions per slice

% returns:
% userDemands: (number of resources x users x time)
% c_uD: (users x time), capacity perceived by each user.

seed=1;
backhaulCapacities = 3;
%% Network settings
users=bsN * sectors * sat; % total number of users
[ NetSettings ] = Network_Settings(sat, bsN, interdistance, users,...
    simulationTime, warmup, backhaulCapacities);
%% Operator settings
% equal shares
belonging= get_belonging(shareVec, users, operators);
shareVec = belonging / users; % readjust the share allocation so that it's
                              % strictly proportional to nUsers.
[ OpSettings ] = Operators_Settings(operators, shareVec, belonging, ...
    NetSettings);
% priorities based on levels
OpSettings.phi_levels = phiLevels;
phi= get_priorities(OpSettings,phiLevels,NetSettings.users);
OpSettings.phi=phi;
OpSettings.alphas=alphas;
%% Scenario
[ bsPositions ] = Scenario_Generation(NetSettings);
%% Mobility
disp('Starting mobility...')
userEquipmentHeight = 1.5;
rad = 500; % m
Mspeed = 1.5; %m/s
[uX, uY, uZ] = RWP_border_circle(NetSettings.users, NetSettings.simulation_time,...
                            userEquipmentHeight, rad, Mspeed);
trace(:,:,1) = uX;
trace(:,:,2) = uY;
disp('done mobility.')
%% Link estimation
disp('Starting Link estimation...')
[ capacityTrace ] = Link_Estimation(NetSettings, trace, bsPositions, 2);
userDemands = zeros(NetSettings.R, NetSettings.users, size(capacityTrace, 3));
for t=1:size(capacityTrace,3)
    [c_u, bs] = max(capacityTrace(:,:,t)');
    c_uD(:, t) = ones(size(c_u)); % to keep it simple, uniform capacity for now.
    %c_uD(:,t) = c_u / 1000000; % in Mbps
    bsD(:,t) = bs; % bs associate
    % also associate the backhaul to form the demanding vectors
    for u = 1:NetSettings.users
        backhaul = bsPositions(ceil(bs(u)/3), 4);
        cDemand = zeros(NetSettings.R, 1);
        cDemand(ceil(bs(u)/3)) = 1;
        cDemand(backhaul) = 1;
        userDemands(:, u, t) = cDemand;
    end
end
OpSettings.c_u=c_uD;
OpSettings.bs=bsD;
disp('done link estimation.')

