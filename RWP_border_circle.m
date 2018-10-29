%% RWP_border_circle
% each milisecond update
function [uX,uY,uZ]=RWP_border_circle(users,dur,UE_height,rad,Mspeed)
duration=dur;
uX=zeros(users,duration);
uY=zeros(users,duration);
uZ=zeros(users,duration);
for u=1:users
    % generate init
    r=sqrt(rand(1,1))*(rad+20);thI=rand()*2*pi;xI=r*[cos(thI)];yI=r*[sin(thI)];
    t=2;
    x=zeros(2e5,1);y=zeros(2e5,1);
    x(1)=xI;y(1)=yI;
    while t<duration
        % generate dest
        r=rad+20;thO=rand()*2*pi;xF=r*[cos(thO)];yF=r*[sin(thO)];
        % pick speed
        sp=Mspeed+1*rand()-0.5;
        % calculate distance in meters
        d=sqrt((xF-xI)^2+(yF-yI)^2);
        % differentials per 1 millisec
        dX=(xF-xI)/floor(d/sp);
        dY=(yF-yI)/floor(d/sp);
        % calculate time to complete and update
        x(t:t+floor(d/sp))=x(t-1)+[dX.*[1:1+floor(d/sp)]]';
        y(t:t+floor(d/sp))=y(t-1)+[dY.*[1:1+floor(d/sp)]]';
        t=t+floor(d/sp);
        xI=xF;yI=yF;
    end
    uX(u,:)=x(1:duration);uY(u,:)=y(1:duration);uZ(u,:)=UE_height*ones(duration,1);
end
end


