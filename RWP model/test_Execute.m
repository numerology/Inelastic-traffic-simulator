function [trace] = rwpmodel(users,xmax,xmin,ymax,ymin)
users=200
s_input = struct('V_POSITION_X_INTERVAL',[xmin xmax],...%(m)
                 'V_POSITION_Y_INTERVAL',[ymin ymax],...%(m)
                 'V_SPEED_INTERVAL',[1.2 10.2],...%(m/s)
                 'V_PAUSE_INTERVAL',[0 1],...%pause time (s)
                 'V_WALK_INTERVAL',[2.00 10.00],...%walk time (s)
                 'V_DIRECTION_INTERVAL',[-180 180],...%(degrees)
                 'SIMULATION_TIME',130,...%(s)
                 'NB_NODES',users);
s_mobility = Generate_Mobility(s_input);
% trace=
%trace contains: (user,time,x/y)
for i=1:users
    trace(i,:,1)=s_mobility.VS_NODE(i).V_POSITION_X;
    trace(i,:,2)=s_mobility.VS_NODE(i).V_POSITION_Y;
end
timeStep = 0.1;%(s)
test_Animate(s_mobility,s_input,timeStep);