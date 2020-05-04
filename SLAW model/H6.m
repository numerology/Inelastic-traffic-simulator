%% H6
% Fixed
num_user=190 * 3;
size_max=1000;
% Hours of simulation time. The simulation time in GREETsimulation is by
% minute, so Thours = simulationTime / 60.
Thours= 5000 / 60; % Max simulation time = 5000
MIN_PAUSE=0;
MAX_PAUSE=30;
beta=1;
% i dependent
n_wp=75;%50+450*i;
B_range=60;%25+250*i;
dist_alpha=4;%1.5+4.5*i;
v_Hurst=0.75;%0.55+0.4*i;
for seed=11:14
        % Generate user movements
        trace=SLAW_MATLAB(dist_alpha, num_user, size_max, n_wp, v_Hurst,...
            Thours, B_range, beta, MIN_PAUSE, MAX_PAUSE);
        save(strcat('./Heterogeneity/H6_long_','seed',num2str(seed)))
end