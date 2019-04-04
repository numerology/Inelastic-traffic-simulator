% mainminreq The script used to simulate the settings where minimal rate
% requirements are specified for users.
% Author: numerology
% 4/4/2019
clc, close all, clear all, gcp;
%% Settings
nSlice = 3;

simulationTime = 10000;
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
threshold = 0.5 * capacity / perBSLoad; % min rate requirement
netSettings = [];
netSettings.bsNS = nBaseStations;
opSettings = [];
opSettings.s_o = shareVec;
opSettings.shareDist = shareDist;
