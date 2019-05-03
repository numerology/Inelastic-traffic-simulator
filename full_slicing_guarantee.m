% Script for SCG simulation
% with more realistic mobility model in addition to poisson.
clc, close all, clear all
%% Set up
nSlices = 4; % num of slices
sat = 3; % U/B (use only integers...)
simulationTime = 1000; % seconds
warmup = 0;
bsN = 19;
sectors = 3;
interdistance = 1000;
% User mobility patterns:
% RWP for roughly uniform spatial loads.
model = {'RWP'}; 
shareVec = 1/3 * ones(1,3); % shares
gcp;

%% Mobility and Link estimation
[NetSettings, OpSettings, capacityPerUser, bs, userPos, bsPos] = ...
    networkconfiguration(simulationTime, ...
    warmup, bsN, sectors,...
    interdistance, model,...
    shareVec, phiLevels, sat, o, alphas, 3);

%% Adjust share distribution for new proposed scheme according to the load distribution
% the sum of share across BSs <= share * |B| per slice.
loadDist = getloaddistribution(OpSettings, NetSettings, bs, simulationTime);
% use a similar heuristic to allocate shares
% OpSettings.shareDist = getsharedistribution(OpSettings, loadDist);
OpSettings.shareDist = loadDist;
OpSettings.s_o = [sum(OpSettings.shareDist, 2)]';