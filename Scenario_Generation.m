%% Scenario Generation
% Pablo Caballero Garcés and JxZheng
% 30/03/15
% Returns:
% (:,1) - X positions
% (:,2) - Y positions
% (:,3) - idx of associated backhaul resource

function [ bsPositions ] = Scenario_Generation(NetSettings)
% Generate front edge base station position,
% together with the association with backhaul resources
% Currently assume BS = 19
%% Generate multiple
assert(NetSettings.bsN == 19)
[bs] = Scenario(NetSettings.interdistance);
bsPositions = bs(1:19, 1:3);
backhaulAssociation = [1 1 1 2 3 3 4 2 2 3 4 4 5 5 7 7 7 6 5]';
bsPositions = [bsPositions backhaulAssociation];
end
