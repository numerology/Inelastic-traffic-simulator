%% Mobility Model Generation
% Pablo Caballero Garcés
% 30/03/15
% (user,time,position)
function [ trace ] = Mobility_Model(model, seed, NetSettings)
maxminD=[-2*173 2*173 -100 300;-2*173 2*173 -300 300;-2*173 2*173 -300 500;
    -3*173 2*173 -300 500;-3*173 2*173 -500 500;-3*173 3*173 -500 500;];
d=[];d(4)=1;d(7)=2;d(10)=3;d(13)=4;d(16)=5;d(19)=6;
Mx=maxminD(d(NetSettings.bsN),2);
mx=maxminD(d(NetSettings.bsN),1);
My=maxminD(d(NetSettings.bsN),4);
my=maxminD(d(NetSettings.bsN),3);

if strcmp(model,'SLAW')
    
    multSeed(:,:)=[Mx/1000 My/1000];
    counter=d(NetSettings.bsN);
    dev(:,:)=[-500,-400;-500,-500;-500,-350;-600,-350;-600,-500;-500,-500];
    valD='NC';dista=2;
    
    % (load precached SLAW model)
    load(strcat('./Models',valD,'/trace_d_',num2str(dista),...
        '_s_',num2str(seed),'.mat'));
    
    trace(:,:,1)=multSeed(:,1)*(trace(:,:,1)+dev(counter,1)); % Adapt range [0,1000] -> [-500,500]
    trace(:,:,2)=multSeed(:,2)*(trace(:,:,2)+dev(counter,2)); % Adapt range [0,1000] -> [-500,500]
elseif strcmp(model,'SLAWH')
    
    multSeed(:,:)=[Mx/900 My/900];
    counter=d(NetSettings.bsN);
    dev(:,:)=[-500,-400;-500,-500;-500,-350;-600,-350;-600,-500;-500,-500];
    valD='NC';dista=2;
    hlevel=NetSettings.hlevel;
    % (load precached SLAW model)
    load(strcat('./SLAW model/Heterogeneity/',hlevel,'/',hlevel,'_seed',...
        num2str(seed),'.mat'));
    
    trace(:,:,1)=multSeed(:,1)*(trace(:,:,1)+dev(counter,1)); % Adapt range [0,1000] -> [-500,500]
    trace(:,:,2)=multSeed(:,2)*(trace(:,:,2)+dev(counter,2)); % Adapt range [0,1000] -> [-500,500]    
elseif strcmp(model,'RWP')
    
    addpath('RWP model/')
    trace=rwpmodel(NetSettings.users,(NetSettings.interdistance/200)*Mx,...
        (NetSettings.interdistance/200)*mx,...
        (NetSettings.interdistance/200)*My,...
        (NetSettings.interdistance/200)*my,...
        NetSettings.warm_up+NetSettings.simulation_time);
    
else
    error('NO MODEL')
end

end