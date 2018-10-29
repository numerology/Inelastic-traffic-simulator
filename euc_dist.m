%% Euclid distance and angle computation

function [ distance,angle ] = euc_dist( x1,y1,x2,y2 )
    distance=sqrt((x2-x1).^2+(y2-y1).^2);
    angle=atan2(x2-x1,y2-y1)*180/pi;
end

