function [NetSettings, OpSettings, capacityPerUser, bsAssociation, trace, ...
    bsPositions]=networkconfiguration(simulationTime, warmup, nBS, nSectors,...
    interdistance, model, shareVec, phiLevels, sat, operators, alphas, ...
    mobilityConfiguration)
%% Network settings
users=nBS * nSectors * sat; % total number of users
[ NetSettings ] = Network_Settings(sat, nBS, interdistance, users, ...
    simulationTime, warmup);
NetSettings.model = model;
%% Operator settings
% equal shares
belonging= get_belonging(shareVec,users,operators);
shareVec=belonging/users;
[ OpSettings ] = Operators_Settings(operators,shareVec,belonging,NetSettings);
% priorities based on levels
OpSettings.phi_levels=phiLevels;
phi= get_priorities(OpSettings,phiLevels,NetSettings.users);
OpSettings.phi=phi;
OpSettings.alphas=alphas;
%% Scenario
[ bsPositions ] = Scenario_Generation(NetSettings);
%% Mobility
disp('Starting mobility...')
UE_height=1.5;rad=100; % m
Mspeed = 10; %m/s

[uX,uY,uZ]=RWP_border_circle(NetSettings.users,NetSettings.simulation_time,...
                             UE_height,rad,Mspeed, 1, OpSettings);
                         
S = load('./SLAW model/Heterogeneity/alpha2nUser1710');
hetTrace = circlewrap(S.trace(1:users, :, 1:2), rad);

S2 = load('./SLAW model/Heterogeneity/H6_seed12');
%S2 = S;
hetTrace2 = circlewrap(S2.trace(1:users, :, 1:2), rad);

S3 = load('./SLAW model/Heterogeneity/H6_seed13');
%S3 = S;
hetTrace3 = circlewrap(S3.trace(1:users, :, 1:2), rad);

S4 = load('./SLAW model/Heterogeneity/H6_seed14');
%S4 = S;
hetTrace4 = circlewrap(S4.trace(1:users, :, 1:2), rad);

nHetTraceUser = size(hetTrace, 1)

switch mobilityConfiguration
    
    case 1 % All RWP
        trace(:, :, 1) = uX(:, :);
        trace(:, :, 2) = uY(:, :);
    case 2 % inelastic SLAW, elastic RWP
        trace(OpSettings.ops_belongs >= 3,:,1) = uX(OpSettings.ops_belongs >= 3, :);
        trace(OpSettings.ops_belongs >= 3,:,2) = uY(OpSettings.ops_belongs >= 3, :);
        for o = 1:2
            trace(OpSettings.ops_belongs == o, 1:NetSettings.simulation_time, :) ... 
                = hetTrace(OpSettings.ops_belongs ...
                == o, 1:NetSettings.simulation_time, :);
        end
    case 3 % All SLAW, with the same hotspots
        trace(:, :, :) = hetTrace(:, 1:NetSettings.simulation_time, :);
    case 4 % All SLAW, with different hotspots
        trace(OpSettings.ops_belongs == 1, 1:NetSettings.simulation_time, :) ... 
            = hetTrace(OpSettings.ops_belongs == 1, 1:NetSettings.simulation_time, :);
        trace(OpSettings.ops_belongs == 2, 1:NetSettings.simulation_time, :) ... 
            = hetTrace(OpSettings.ops_belongs == 2, 1:NetSettings.simulation_time, :);
        trace(OpSettings.ops_belongs == 3, 1:NetSettings.simulation_time, :) ... 
            = hetTrace(OpSettings.ops_belongs == 3, 1:NetSettings.simulation_time, :);
        trace(OpSettings.ops_belongs == 4, 1:NetSettings.simulation_time, :) ... 
            = hetTrace(OpSettings.ops_belongs == 4, 1:NetSettings.simulation_time, :);
        
    case 5 % super heteroneous case, all SLAW with different hotspots, but also 
           % inelastic users are of different phi within a single slice. 
           % its mobility pattern is the same as Scenario 4.
        trace(OpSettings.ops_belongs == 1, 1:NetSettings.simulation_time, :) ... 
            = hetTrace(OpSettings.ops_belongs == 1, 1:NetSettings.simulation_time, :);
        trace(OpSettings.ops_belongs == 2, 1:NetSettings.simulation_time, :) ... 
            = hetTrace(OpSettings.ops_belongs == 2, 1:NetSettings.simulation_time, :);
        trace(OpSettings.ops_belongs == 3, 1:NetSettings.simulation_time, :) ... 
            = hetTrace(OpSettings.ops_belongs == 3, 1:NetSettings.simulation_time, :);
        trace(OpSettings.ops_belongs == 4, 1:NetSettings.simulation_time, :) ... 
            = hetTrace(OpSettings.ops_belongs == 4, 1:NetSettings.simulation_time, :);
end

         


disp('done mobility.')

%% Mobility heterogeneous
% disp('Starting mobility...')
% trace = Mobility_Model('SLAWH', 1, NetSettings);
% disp('done mobility.')
%% Link estimation
disp('Starting Link estimation...')
[ c_ijt ] = LinkEstimation(NetSettings,trace,bsPositions,2);
for t=1:size(c_ijt,3)
    [c_u, bs]=max(c_ijt(:,:,t)');
    capacityPerUser(:,t)=c_u/1000000; % in Mbps
    bsAssociation(:,t)=bs;
end
OpSettings.c_u=capacityPerUser;
OpSettings.bs=bsAssociation;
disp('done link estimation.')

