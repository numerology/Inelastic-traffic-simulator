function [NetSettings, OpSettings, capacityPerUser, bsAssociation, trace, ...
    bsPositions]=networkconfiguration(simulationTime, warmup, nBS, nSectors,...
    interdistance, model, shareVec, phiLevels, sat, operators, alphas, seed)
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
UE_height=1.5;rad=500; % m
Mspeed=1.5; %m/s
[uX,uY,uZ]=RWP_border_circle(NetSettings.users,NetSettings.simulation_time,...
                            UE_height,rad,Mspeed);
trace(:,:,1)=uX;
trace(:,:,2)=uY;
disp('done mobility.')
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

