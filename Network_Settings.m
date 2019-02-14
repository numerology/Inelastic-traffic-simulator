%% Network configuration settings
% Pablo Caballero Garcia
% 30/03/15
function [ NetSettings ] = Network_Settings(saturation, bsN, interdistance, ...
    users, simulationTime, warmup)

    NetSettings=[];
    
    NetSettings.saturation=saturation; % Saturation
    
    NetSettings.bsN=bsN; % base stations

    NetSettings.bsNS=bsN*3; % |B|

    NetSettings.interdistance=interdistance; % Inter site distance
   
    NetSettings.users=users; % |U|
    
    NetSettings.simulation_time=simulationTime; % simulation_time
    
    NetSettings.fc=2.5; % frecuency carrier (Ghz)
    
    NetSettings.warm_up=warmup; % Specify the time of link estimation last time to pick up. (obsolete)

    NetSettings.m=5; % m for Online
    
    NetSettings.diffR=0; % THR for Online
    
end   
