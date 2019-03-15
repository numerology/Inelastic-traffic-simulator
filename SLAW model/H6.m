%% H6
% Fixed
num_user=190 * 3;
size_max=1000;
Thours= 1000 / 60;
MIN_PAUSE=0;
MAX_PAUSE=30;
beta=1;
for seed=11:20
        %% i dependent
        n_wp=75%50+450*i;
        B_range=60%25+250*i;
        dist_alpha=4%1.5+4.5*i;
        v_Hurst=0.75%0.55+0.4*i;
        %% Generate user movements
        trace=SLAW_MATLAB(dist_alpha, num_user, size_max, n_wp, v_Hurst,...
            Thours, B_range, beta, MIN_PAUSE, MAX_PAUSE);
        filename = strcat('./Heterogeneity/H6_','seed',num2str(seed));
        save(filename)
end