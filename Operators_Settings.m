%% Operators configuration settings
% Pablo Caballero Garcés
% 30/03/15
function [ OpSettings ] = Operators_Settings(operators, shareVec, ...
    belonging, NetSettings)
% Generate the profile of network slices/operators
% operators: number of operators, |V|
% shareVec: vector of share allocation among slices
% belonging: vector, number of users under each slice
% NetSettings: network profile.

%% Number of operators
OpSettings.operators = operators; 
assert(size(shareVec, 2) == operators, 'Shares wrong input');
assert(size(belonging, 2) == operators, 'Belonging wrong input');

%% Network shares
OpSettings.shareVec = shareVec; 
%% Users per operator and ops_belong matrix
opsBelongs = [];

for i = 1:operators
    opsBelongs=[opsBelongs; i * ones(belonging(i),1)];
end
OpSettings.N_k = belonging; % vector of number of users per slice.
OpSettings.opsBelongs = opsBelongs'; % user belonging vector, 
% ops_belongs(i) = the idx of slice user i belongs to.

%% Weights
% According to equal weight sharing in SCPF,
for i=1:NetSettings.users
    userWeights(i) = shareVec(opsBelongs(i)) / belonging(opsBelongs(i));
end
OpSettings.w_i = userWeights;    
%% Checks
assert(floor(100*sum(OpSettings.w_i)+0.000001) == 100, 'Weight error');
end
