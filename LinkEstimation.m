%% Link Estimation
% Pablo Caballero Garces, Numerology
% 30/03/15
function [ capacities ] = LinkEstimation(NetSettings, trace, bsPositions, multFading)
% Generate the link capacity given the position and fading
% Parameters:
% Returns:

TXp=41;ThermalNoisePlusInterfig=-104;
capacities = zeros(NetSettings.users, NetSettings.bsNS, ...
    NetSettings.simulation_time);
ppm = ParforProgMon('Link estimating : ', NetSettings.simulation_time);
parfor t = 1:NetSettings.simulation_time
    %if rem(t,200)==0; disp(t);end
    distA=zeros(size(unique(bsPositions(:,3)),2),1);
    angleA=zeros(size(unique(bsPositions(:,3)),2),1);
    tmpCapacities = zeros(NetSettings.users, NetSettings.bsNS);
    for l = 1:NetSettings.users
        %% Distance and angle to each bs
        [dist angle]=euc_dist(trace(l,t,1),trace(l,t,2),bsPositions(:,1),bsPositions(:,2));
        c=0;
        
        for bs=[unique(bsPositions(:,3))']
            c=c+1;
            eqs=find(bsPositions(:,3)==bs);
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
        bsNS = zeros(NetSettings.bsN, 3);
        for bs=1:NetSettings.bsN
            [ gain_sectorized ] = antenna_gain(angle(bs)+180);
            bsNS(bs,:)=gain_sectorized;
        end
        %% Tx power
        P_tx=10^(TXp);
        %% Noise
        alp=db2pow(ThermalNoisePlusInterfig);
        %% Fading
        sh_fad=db2pow(normrnd(0,multFading));
        %% Channel Gain
        g=db2pow(pathlossV+bsNS+17-3);
        g=reshape(g',[NetSettings.bsNS 1]);
        %% Balance
        SI = zeros(1, NetSettings.bsNS);
        for j=1:NetSettings.bsNS
            Srec=(P_tx*g(j));
            Inter=sum(P_tx*g)-(P_tx*g(j));
            SI(j)=Srec*sh_fad/(alp+1*Inter);
        end
        %SINR_db=pow2db(SI);
        c=real(log2(1+SI));
        cqi=c2CQI(c);
        tbs=cqi2tbs(cqi);
        tmpCapacities(l,:) = tbs;
        if tbs==0
            error('SINR error')
        end
    end
    capacities(:,:,t) = tmpCapacities;
    ppm.increment();
end
capacities(capacities==0)=-Inf;
end

