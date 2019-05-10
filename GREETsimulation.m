function [outputArg1,outputArg2] = GREETsimulation(satVector, profile)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
outputArg1 = inputArg1;
outputArg2 = inputArg2;


nSlices = 4; % num of slices

sat = 3; % U/B (use only integers...)
simulationTime = 2000; % seconds

phiLevels = 1;alphas = [1, 1, 1, 1]; % legacy parameters
warmup = 0;
bsN = 19;
sectors = 3;
interdistance = 15;
outageTol = 0.01;
minSharePerBS = 0.001;
% User mobility patterns:
% RWP for roughly uniform spatial loads.
model = {'RWP'};
shareVec = 1/4 * ones(1, 4); % this only means user numbers are the same.
sliceCats = [0 0 1 1];


end