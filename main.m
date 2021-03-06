clc, close all, clear all
%% Settings
o = 3; % num of slices
sat = 1; % U/B (use only integers...)
simulationTime = 1000; % seconds
phiLevels = 1;alphas = [1, 1, 1];
warmup = 0;bsN = 19;sectors = 3;
interdistance = 1000;
% User mobility patterns:
% RWP for roughly uniform spatial loads.
model = {'RWP'}; 
shareVec = 1/3 * ones(1,3); % shares
gcp;

%% Mobility and Link estimation
[NetSettings, OpSettings, capacityPerUser, bs, userPos, bsPos] = ...
    networkconfiguration(simulationTime, ...
    warmup, bsN, sectors,...
    interdistance, model,...
    shareVec, phiLevels, sat, o, alphas, 3);

%% Adjust share distribution for new proposed scheme according to the load distribution
% the sum of share across BSs <= share * |B| per slice.
loadDist = getloaddistribution(OpSettings, NetSettings, bs, simulationTime);
% use a similar heuristic to allocate shares
% OpSettings.shareDist = getsharedistribution(OpSettings, loadDist);
OpSettings.shareDist = loadDist;
OpSettings.s_o = [sum(OpSettings.shareDist, 2)]';
%% Compute fractions
ppm = ParforProgMon('Simulating resource sharing : ', NetSettings.simulation_time);
parfor t=1:simulationTime
   
%     [r,f,b] = Static_Slicing(NetSettings, OpSettings, [capacityPerUser(:,t)]', [bs(:,t)]');
%     rates_SS(:,t)=r;
%     fractions_SS(:,t)=f;
%     btd_SS(:,t)=b;
%     [r,f,b] = flexibleSCPF(NetSettings, OpSettings, [capacityPerUser(:,t)]', [bs(:,t)]');
%     rates_fSCPF(:,t)=r;
%     fractions_fSCPF(:,t)=f;
%     btd_fSCPF(:,t)=b;
    [r,f,b] = Static_Slicing(NetSettings, OpSettings, [capacityPerUser(:,t)]', [bs(:,t)]');
    rates_GPS(:,t)=r;
    fractions_GPS(:,t)=f;
    btd_GPS(:,t)=b;
    [r,f,b] = SCPF(NetSettings, OpSettings, [capacityPerUser(:,t)]', [bs(:,t)]');
    rates_SCPF(:,t)=r;
    fractions_SCPF(:,t)=f;
    btd_SCPF(:,t)=b;
    [r,f,b] = DIFFPRICE(NetSettings, OpSettings, [capacityPerUser(:,t)]', [bs(:,t)]');
    rates_DIFFPRICE(:,t)=r;
    fractions_DIFFPRICE(:,t)=f;
    btd_DIFFPRICE(:,t)=b;
    [r,f,b] = MAXWEIGHT(NetSettings, OpSettings, [capacityPerUser(:,t)]', [bs(:,t)]');
    rates_MAXWEIGHT(:,t)=r;
    fractions_MAXWEIGHT(:,t)=f;
    btd_MAXWEIGHT(:,t)=b;
    ppm.increment();
end
%% Plot performance seen by some randomly selected users.
i1=30;
i2=113;
i3=243;
figure();
subplot(3,1,1)
plot(btd_fSCPF(i1,:),'-.b')
hold on
plot(btd_GPS(i1,:),'--r')
plot(btd_SCPF(i1,:),'-g')
plot(btd_capSCG(i1,:),':k')
subplot(3,1,2)
plot(btd_fSCPF(i2,:),'-.b')
hold on
plot(btd_GPS(i2,:),'--r')
plot(btd_SCPF(i2,:),'-g')
plot(btd_capSCG(i2,:),':k')
subplot(3,1,3)
plot(btd_fSCPF(i3,:),'-.b')
hold on
plot(btd_GPS(i3,:),'--r')
plot(btd_SCPF(i3,:),'-g')
plot(btd_capSCG(i3,:),':k')
legend('SCG', 'GPS','SCPF','SCG')
%% Take a look at the mean performance
disp('Overall')
fprintf('mean btd of GPS = %f\n', mean(mean(btd_GPS)));
fprintf('mean btd of SCPF = %f\n', mean(mean(btd_SCPF)));
fprintf('mean btd of DIFFPRICE = %f\n', mean(mean(btd_DIFFPRICE)));
fprintf('mean btd of MAXWEIGHT = %f\n', mean(mean(btd_MAXWEIGHT)));
%% Take a look at the mean performance for a specific slice
for sliceIdx = 1:o
    fprintf('For slice %i\n', sliceIdx);
    fprintf('mean btd of GPS = %f\n', ...
        mean(mean(btd_GPS(OpSettings.ops_belongs == sliceIdx, :, :))));
    fprintf('mean btd of SCPF = %f\n', ...
        mean(mean(btd_SCPF(OpSettings.ops_belongs == sliceIdx, :, :))));
    fprintf('mean btd of MAXWEIGHT = %f\n', ...
        mean(mean(btd_MAXWEIGHT(OpSettings.ops_belongs == sliceIdx, :, :))));
    fprintf('mean btd of DIFFPRICE = %f\n', ...
        mean(mean(btd_DIFFPRICE(OpSettings.ops_belongs == sliceIdx, :, :))));
end

% fprintf('mean rate of SCPF = %f\n', ...
%     mean(mean(rates_SCPF(OpSettings.ops_belongs == sliceIdx, :, :))));
% fprintf('mean rate of SCG = %f\n', ...
%     mean(mean(rates_capSCG(OpSettings.ops_belongs == sliceIdx, :, :))));
% fprintf('mean rate of GPS = %f\n', ...
%     mean(mean(rates_GPS(OpSettings.ops_belongs == sliceIdx, :, :))));
% fprintf('mean rate of Flexible SCPF = %f\n', ...
%     mean(mean(rates_fSCPF(OpSettings.ops_belongs == sliceIdx, :, :))));
%% Some CDF of BTD plot
figure()
cdfplot(reshape(log(rates_GPS), [1, size(btd_GPS, 1) * size(btd_GPS, 2)]));
hold on
cdfplot(reshape(log(rates_SCPF), [1, size(btd_SCPF, 1) * size(btd_SCPF, 2)]));
cdfplot(reshape(log(rates_MAXWEIGHT), [1, size(btd_SCPF, 1) * size(btd_SCPF, 2)]));
cdfplot(reshape(log(rates_DIFFPRICE), [1, size(btd_SCPF, 1) * size(btd_SCPF, 2)]));
title('CDF of log BTD')
xlabel('BTD')
%ylim([0.9 1]);
legend('SS', 'SCPF', 'MAXWEIGHT', 'DIFFPRICE');
%% 
pl=0;
if pl==1
%% Some plotting all users at a given time
Scenario_draw(200);
plot(bsPos(1:19,1),bsPos(1:19,2),'k^')
hold on
plot(userPos(:,1000,1),userPos(:,1000,2),'rs')
axis equal
xlim([-500, 500])
ylim([-500, 500])
hold off
%% Some plotting for 5 users
figure();
Scenario_draw(200);
plot(bsPos(1:19,1),bsPos(1:19,2),'k^')
axis equal
hold on
for t=1:simulationTime
    t
    plot(userPos(1,t,1),userPos(1,t,2),'s','Color',[1/4+3/4*(capacityPerUser(6,t)/max(capacityPerUser(:))),0,0])
    plot(userPos(2,t,1),userPos(2,t,2),'s','Color',[0,1/4+3/4*(capacityPerUser(2,t)/max(capacityPerUser(:))),0])
    plot(userPos(3,t,1),userPos(3,t,2),'s','Color',[0,0,1/4+3/4*(capacityPerUser(3,t)/max(capacityPerUser(:)))])
    plot(userPos(4,t,1),userPos(4,t,2),'s','Color',[1/4+3/4*(capacityPerUser(4,t)/max(capacityPerUser(:)))...
                                                       ,1/4+3/4*(capacityPerUser(4,t)/max(capacityPerUser(:))),0])
    plot(userPos(5,t,1),userPos(5,t,2),'s','Color',[1/4+3/4*(capacityPerUser(5,t)/max(capacityPerUser(:)))...
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