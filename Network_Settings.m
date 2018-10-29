%% Network configuration settings
% Pablo Caballero Garcés
% 30/03/15
function [ NetSettings ] = Network_Settings(saturation,bsN,interdistance,users,simulationTime,warmup, backhaulCapacities)

    NetSettings=[];
    
    NetSettings.saturation=saturation; % Saturation
    
    NetSettings.bsN=bsN; % base stations |B|, only including the front
    
    NetSettings.R = bsN + 7; % currently only study the 19 front with 7 backhaul setting.

    NetSettings.interdistance=interdistance; % Inter site distance
   
    NetSettings.users=users; % |U|
    
    NetSettings.simulation_time=simulationTime; % simulation_time
    
    NetSettings.fc=2.5; % frecuency carrier (Ghz)
    
    NetSettings.warm_up=warmup; % frecuency carrier (Ghz)

    NetSettings.m=5; % m for Online
    
    NetSettings.diffR=0; % THR for Online
    
    NetSettings.backhaulCapacity = backhaulCapacities;
    
end   
