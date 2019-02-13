clc, close all, clear all
%% Settings
o=3; % num of slices
sat=30; % U/B (use only integers...)
simulationTime=5; % seconds
shareVec=[1/3 1/3 1/3]; % shares
%%
phiLevels=1;alphas=[1,1,1];
warmup=0;bsN=19;sectors=3;
interdistance=200; t=10; model={'RWP'};
%% Mobility and Lik estimation
[NetSettings, OpSettings, capacityPerUser, bs, usersPos,bsPosition, ...
    userDemands]=Network_configuration(simulationTime, warmup,bsN,sectors,...
    interdistance, shareVec, phiLevels, sat, o, alphas, 0);

%% Compute fractions
for t=1:simulationTime
    t
    [rate,fraction,btd]=StaticSlicing(NetSettings, OpSettings, ...
        capacityPerUser(:,t)' ,bs(:,t)', userDemands(:, :, t));
    ratesSS(:,t)=rate;
    fractionsSS(:,t)=fraction;
    btdSS(:,t)=btd;
%     [rate,fraction,btd]=GPS(NetSettings, OpSettings,capacityPerUser(:,t)',bs(:,t)');
%     ratesGPS(:,t)=rate;
%     fractionsGPS(:,t)=fraction;
%     btdGPS(:,t)=btd;
%     [rate,fraction,btd]=SCPF(NetSettings, OpSettings,capacityPerUser(:,t)',bs(:,t)');
%     ratesSCPF(:,t)=rate;
%     fractionsSCPF(:,t)=fraction;
%     btdSCPF(:,t)=btd;
end
%%
i1=1;
i2=2;
i3=3;
subplot(3,1,1)
plot(ratesSS(i1,:),'-.')
hold on
plot(ratesGPS(i1,:),'--')
plot(ratesSCPF(i1,:),'-')
subplot(3,1,2)
plot(ratesSS(i2,:),'-.')
hold on
plot(ratesGPS(i2,:),'--')
plot(ratesSCPF(i2,:),'-')
subplot(3,1,3)
plot(ratesSS(i3,:),'-.')
hold on
plot(ratesGPS(i3,:),'--')
plot(ratesSCPF(i3,:),'-')
legend('SS','GPS','SCPF')
%%
pl=0;
if pl==1
%% Some plotting all users at a given time
Scenario_draw(200);
plot(bsPosition(1:19,1),bsPosition(1:19,2),'k^')
hold on
plot(usersPos(:,1000,1),usersPos(:,1000,2),'rs')
axis equal
xlim([-500, 500])
ylim([-500, 500])
hold off
%% Some plotting for 5 users
Scenario_draw(200);
plot(bsPosition(1:19,1),bsPosition(1:19,2),'k^')
axis equal
hold on
for t=1:simulationTime
    t
    plot(usersPos(1,t,1),usersPos(1,t,2),'s','Color',[1/4+3/4*(capacityPerUser(6,t)/max(capacityPerUser(:))),0,0])
    plot(usersPos(2,t,1),usersPos(2,t,2),'s','Color',[0,1/4+3/4*(capacityPerUser(2,t)/max(capacityPerUser(:))),0])
    plot(usersPos(3,t,1),usersPos(3,t,2),'s','Color',[0,0,1/4+3/4*(capacityPerUser(3,t)/max(capacityPerUser(:)))])
    plot(usersPos(4,t,1),usersPos(4,t,2),'s','Color',[1/4+3/4*(capacityPerUser(4,t)/max(capacityPerUser(:)))...
                                                       ,1/4+3/4*(capacityPerUser(4,t)/max(capacityPerUser(:))),0])
    plot(usersPos(5,t,1),usersPos(5,t,2),'s','Color',[1/4+3/4*(capacityPerUser(5,t)/max(capacityPerUser(:)))...
                                                     ,0,1/4+3/4*(capacityPerUser(5,t)/max(capacityPerUser(:)))])
    xlim([-540, 540])
    ylim([-540, 540])
    pause(eps)
end
disp('done')
end


% 
% 
% 
% 
% 
% 
% %% Utility functions
% for o=1:OpSettings.operators
%     al=OpSettings.alphas(o);
%     if al==1
%         Ut_o{o}=@(x) sum(OpSettings.phi(OpSettings.ops_belongs==o).*log(x));
%     else
%         Ut_o{o}=@(x) sum(OpSettings.phi(OpSettings.ops_belongs==o).*(x).^(1-al)./(1-al));
%     end
% end
% %% Define initial weight state (random ones)
% w0=[];
% for i=1:OpSettings.operators
%     wi=rand(1,OpSettings.N_k(i));swi=sum(wi,2);
%     w0=[w0 ,OpSettings.s_o(i)*wi(1,:)/swi];
% end
% %% Global
% %
% disp('Starting global...')
% if algorithms(1)==1,[weights_M,U_M]=Social_optimum(NetSettings, OpSettings,c_u,bs,extra);end
% U_M
% %% SS
% disp('Starting SS...')
% if algorithms(2)==1,[weights_S,U_S,U_S_o]=Static_Slicing(NetSettings, OpSettings,c_u,bs,extra,Ut_o);end
% U_S
% %% Nash EQ
% disp('Starting Dist...')
% preround=0;phi=OpSettings.phi;weights=phi;it=0;
% round=0;
% while (sum(abs(preround-weights))>3e-12)
%     round=round+1;
%     disp('____________________________________________________')
%     preround=weights;
%     for itd=1:OpSettings.operators
%         it=it+1;w_pre=weights;
%         o=1+rem(it-1,OpSettings.operators);
%         alpha=OpSettings.alphas(o);
%         if alpha>=1
%             trig=1;
%             
%             while trig>1e-14
%                 m=weights(OpSettings.ops_belongs==o);
%                 b_m=bs(OpSettings.ops_belongs==o);
%                 p_m=phi(OpSettings.ops_belongs==o);
%                 down=zeros(size(m,2),1);
%                 for u=1:size(m,2)
%                     b_m(u);
%                     a=sum(weights(OpSettings.ops_belongs~=o & bs==b_m(u)));
%                     l=sum(weights(OpSettings.ops_belongs==o & bs==b_m(u)));
%                     if a+l>0
%                         down(u)=c_u(u)^(1/alpha-1)*p_m(u)^(1/alpha)*(a)^(1/alpha)/(a+l)^(2/alpha-1);
%                     end
%                 end
%                 for u=1:size(m,2)
%                     a=sum(weights(OpSettings.ops_belongs~=o & bs==b_m(u)));
%                     l=sum(weights(OpSettings.ops_belongs==o & bs==b_m(u)));
%                     m_post(u)=OpSettings.s_o(o)*c_u(u)^(1/alpha-1)*p_m(u)^(1/alpha)*(a^(1/alpha)/(a+l)^(2/alpha-1))/sum(down);
%                 end
%                 m_post(isnan(m_post))=0;
%                 m_post(m_post==0)=eps;
%                 weights(OpSettings.ops_belongs==o)=m_post;
%                 trig=sum(abs(m_post-m));
%             end
%             clear m m_post
%         else
%             m=weights(OpSettings.ops_belongs==o);
%             [ U_1 ] = ut_op(m,OpSettings,NetSettings,alpha,c_u,weights,o,bs);
%             opts1=optimoptions('fmincon');
%             opts1.Display = 'off';
%             %opts1.MaxFunEvals = 5e4;
%             opts1.TolFun=min(abs(U_M*1e-2),1e-4);
%             [x U_2] = fmincon(@(x)-ut_op(x,OpSettings,NetSettings,alpha,c_u,weights,o,bs),...
%                 m,...
%                 [],[],...
%                 ones(1,OpSettings.N_k(o)),OpSettings.s_o(o),...
%                 zeros(1,OpSettings.N_k(o)),ones(OpSettings.N_k(o),1),[],opts1);
%             pre=weights;
%             weights(OpSettings.ops_belongs==o)=x;
%             U_T_1=ut_global(weights,OpSettings,NetSettings,alpha,c_u,bs,1);
%             %disp(strcat(num2str(U_1),'<=',num2str(-U_2),'- -',num2str(U_1<-U_2),'- - :',num2str(sum(abs(pre-weights)))))
%         end
%         
%     end
%     U_T_1=ut_global(weights,OpSettings,NetSettings,alpha,c_u,bs,1);
%     disp(strcat(num2str(U_S),'<=', num2str(U_T_1),'<=',...
%         num2str(U_M),'- ',num2str(U_T_1>U_M),...
%         '- ',num2str(U_T_1>U_S),...
%         '- ',num2str(sum(abs(preround-weights)))...
%         ));
%     weights_D=weights;U_D=U_T_1;
% end
% round
% bar([U_S,U_D,U_M])
% 
% if extra_compute==1
% if algorithms(2)==1
%     disp('Starting extra SS...'),extra=[1 2];
%     for i=1:10
%         i
%         e=mean(extra);
%         
%         [~,U_S_e,~]=Static_Slicing(NetSettings, OpSettings,c_u,bs,e);
%         
%         if U_S_e>U_D
%             extra(2)=e;
%         else
%             extra(1)=e;
%         end
%     end
%     extra_S=mean(extra)
% end
% 
% if algorithms(1)==1
%     disp('Starting extra Dist...'),extra=[1 2];
%     for i=1:10
%         e=mean(extra);
%         q=ut_global(weights_D,OpSettings,NetSettings,alpha,c_u,bs,e);
%         if q>U_M
%             extra(2)=e;
%         else
%             extra(1)=e;
%         end
%     end
%     extra_D=mean(extra)
% end
% end
% if ef_compute==1
% for b=1:NetSettings.bsNS
% price(b)=sum(weights_D(bs==b));
% end
% for o=1:OpSettings.operators
%     for u=1:OpSettings.N_k(1)
%         w_tilde(u)=(OpSettings.phi(u)*sum(weights_D(bs==bs(u)&OpSettings.ops_belongs==o)))...
%             /sum(phi(bs==bs(u)&OpSettings.ops_belongs==1));
%         if w_tilde(u)==0
%             w_tilde(u)=eps;
%         end
%         r_tilde(u)=c_u(u)*w_tilde(u)/price(bs(u));
%     end
%     pause()
%     uu(o)= Ut_o{o}(r_tilde);
% end
% ef=max(uu(2:end))-uu(1);
% envy_free=(ef);
% end
% toc
% save(name)