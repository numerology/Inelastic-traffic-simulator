function [trace] = rwpmodel(users,xmax,xmin,ymax,ymin,time)
s_input = struct('V_POSITION_X_INTERVAL',[xmin xmax],...%(m)
    'V_POSITION_Y_INTERVAL',[ymin ymax],...%(m)
    'V_SPEED_INTERVAL',[1.2 10.2],...%(m/s)
    'V_PAUSE_INTERVAL',[0 1],...%pause time (s)
    'V_WALK_INTERVAL',[2.00 40.00],...%walk time (s)
    'V_DIRECTION_INTERVAL',[-180 180],...%(degrees)
    'SIMULATION_TIME',time,...%(s)
    'NB_NODES',users);
s_mobility = Generate_Mobility(s_input);
trace=zeros(users,time,3);
for i=1:users
    for t=1:time
        m=max(find(round(s_mobility.VS_NODE(i).V_TIME)<t));
        trace(i,t,1)=s_mobility.VS_NODE(i).V_POSITION_X(m);
        trace(i,t,2)=s_mobility.VS_NODE(i).V_POSITION_Y(m);
        trace(i,t,4)=sqrt(s_mobility.VS_NODE(i).V_SPEED_X(m)+s_mobility.VS_NODE(i).V_SPEED_Y(m));
    end
end
% timeStep = 1;%(s)
% test_Animate(s_mobility,s_input,timeStep);