% script to get the tradeoff comparison among the 4 benchmarks.
clc, close all, clear all
parpool('local', 40);
warning('off','all');

satVec = [1 3 5 7 10 15 20]; % subject to what type of SLAW model to use
simulationTime = 1000;

poutageGain = zeros(4, length(satVec));
utilGain = zeros(4, length(satVec));

for netProfile = 1:4
    fprintf('beginning configuration no. %d\n', netProfile);
    [poutageGain(netProfile, :), utilGain(netProfile, :)] = ...
        GREETsimulation(satVec, simulationTime, netProfile);
end