%% Scenario Generation
% Pablo Caballero Garcés
% 30/03/15
% Returns:
% (:,1) - X positions
% (:,2) - Y positions
% (:,3) - Mirror equivalence

function [ bs_positions ] = Scenario_Generation(NetSettings)

%% Generate multiple
[bs] = Scenario(NetSettings.interdistance);

%% Mirror system

if NetSettings.bsN==19
    bs(:,3)=bs(:,3);
elseif NetSettings.bsN==16
    bs(bs(:,3)==19,3)=15; bs(bs(:,3)==14,3)=10; bs(bs(:,3)==9,3)=5;
elseif NetSettings.bsN==13
    bs(bs(:,3)==19,3)=15; bs(bs(:,3)==14,3)=10; bs(bs(:,3)==9,3)=5;
    bs(bs(:,3)==16,3)=2; bs(bs(:,3)==17,3)=1; bs(bs(:,3)==18,3)=4;
elseif NetSettings.bsN==10
    bs(bs(:,3)==19,3)=12; bs(bs(:,3)==14,3)=11; bs(bs(:,3)==9,3)=6;
    bs(bs(:,3)==16,3)=2; bs(bs(:,3)==17,3)=1; bs(bs(:,3)==18,3)=4;
    bs(bs(:,3)==15,3)=13; bs(bs(:,3)==10,3)=8; bs(bs(:,3)==5,3)=3;
elseif NetSettings.bsN==7
    bs(bs(:,3)==19,3)=12; bs(bs(:,3)==14,3)=11; bs(bs(:,3)==9,3)=6;
    bs(bs(:,3)==16,3)=6; bs(bs(:,3)==17,3)=3; bs(bs(:,3)==18,3)=8;
    bs(bs(:,3)==15,3)=13; bs(bs(:,3)==10,3)=8; bs(bs(:,3)==5,3)=3;
    bs(bs(:,3)==1,3)=12; bs(bs(:,3)==2,3)=11; bs(bs(:,3)==4,3)=13;
elseif NetSettings.bsN==4
    bs(bs(:,3)==19,3)=7; bs(bs(:,3)==14,3)=6; bs(bs(:,3)==9,3)=6;
    bs(bs(:,3)==16,3)=6; bs(bs(:,3)==17,3)=3; bs(bs(:,3)==18,3)=8;
    bs(bs(:,3)==15,3)=8; bs(bs(:,3)==10,3)=8; bs(bs(:,3)==5,3)=3;
    bs(bs(:,3)==1,3)=3; bs(bs(:,3)==2,3)=6; bs(bs(:,3)==4,3)=8;
    bs(bs(:,3)==11,3)=6; bs(bs(:,3)==12,3)=7; bs(bs(:,3)==13,3)=7;
else
    error('Non programmed mirror for this number of bs')
end

%% Compile
bs_positions=bs;

end
