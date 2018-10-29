%% Link Estimation
% Pablo Caballero Garces
% 30/03/15
function [ c_ijt ] = Link_Estimation(NetSettings,trace,bs_positions,mult_fading)
t=0;
TXp=41;ThermalNoisePlusInterfig=-104;
for tictime=NetSettings.warm_up+1:1:NetSettings.warm_up+NetSettings.simulation_time
    t=t+1;
    fprintf('.');
    if rem(t,100)==0;fprintf(repmat('\b', 1, 100));
        fprintf('%d %%\n',round(100*t/NetSettings.simulation_time));
    end
    %if rem(t,200)==0; disp(t);end
    distA=zeros(size(unique(bs_positions(:,3)),2),1);
    angleA=zeros(size(unique(bs_positions(:,3)),2),1);
    for l=1:NetSettings.users
        %% Distance and angle to each bs
        [dist angle]=euc_dist(trace(l,tictime,1),trace(l,tictime,2),bs_positions(:,1),bs_positions(:,2));
        c=0;
        
        for bs=[unique(bs_positions(:,3))']
            c=c+1;
            eqs=find(bs_positions(:,3)==bs);
            [distA(c), index]=min(dist(eqs));
            angleA(c)=angle(eqs(index));
        end
        dist=distA; angle=angleA;
        %% Path loss
        pathlossV=zeros(NetSettings.bsN,3);
        %angle_Dop=[60-(angle+180); 120-(angle+180); 180-(angle+180)]';
        for bs=1:NetSettings.bsN
            pathlossV(bs,1:3)=-pathloss(dist(bs),NetSettings.fc);
        end
        %% Angle loss
        for bs=1:NetSettings.bsN
            [ gain_sectorized ] = antenna_gain(angle(bs)+180);
            bsNS(bs,:)=gain_sectorized;
        end
        %% Tx power
        P_tx=10^(TXp);
        %% Noise
        alp=db2pow(ThermalNoisePlusInterfig);
        %% Fading
        sh_fad=db2pow(normrnd(0,mult_fading));
        %% Channel Gain
        g=db2pow(pathlossV+bsNS+17-3);
        g=reshape(g',[NetSettings.bsNS 1]);
        %% Balance
        for j=1:NetSettings.bsNS
            Srec=(P_tx*g(j));
            Inter=sum(P_tx*g)-(P_tx*g(j));
            SI(j)=Srec*sh_fad/(alp+1*Inter);
        end
        %SINR_db=pow2db(SI);
        c=real(log2(1+SI));
        cqi=c2CQI(c);
        tbs=cqi2tbs(cqi);
        c_ijt(l,:,t)=tbs;
        if max(c_ijt(l,:,t))==0
            error('SINR error')
        end
    end
end
c_ijt(c_ijt==0)=-Inf;
end

