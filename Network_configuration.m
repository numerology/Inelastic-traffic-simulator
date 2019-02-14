function [NetSettings, OpSettings, c_uD, bsD,trace,bs_positions]=Network_configuration(simulationTime,warmup,bsN,sectors,...
    interdistance, model, s_o, phi_levels, sat, operators, alphas, t)
seed=1;
%% Network settings
users=bsN*sectors*sat; % total number of users
[ NetSettings ] = Network_Settings(sat, bsN, interdistance, users, ...
    simulationTime, warmup);
NetSettings.model = model;
%% Operator settings
% equal shares
belonging= get_belonging(s_o,users,operators);
s_o=belonging/users;
[ OpSettings ] = Operators_Settings(operators,s_o,belonging,NetSettings);
% priorities based on levels
OpSettings.phi_levels=phi_levels;
phi= get_priorities(OpSettings,phi_levels,NetSettings.users);
OpSettings.phi=phi;
OpSettings.alphas=alphas;
%% Scenario
[ bs_positions ] = Scenario_Generation(NetSettings);
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
[ c_ijt ] = LinkEstimation(NetSettings,trace,bs_positions,2);
for t=1:size(c_ijt,3)
    [c_u, bs]=max(c_ijt(:,:,t)');
    c_uD(:,t)=c_u/1000000; % in Mbps
    bsD(:,t)=bs;
end
OpSettings.c_u=c_uD;
OpSettings.bs=bsD;
disp('done link estimation.')

