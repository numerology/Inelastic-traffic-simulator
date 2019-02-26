% Generate clustered user movement trace
%% Users settings
dist_alpha=6;
% need to support at most 19 * 30
num_user = 570 * 3;
size_max = 1000;
simulationTime = 5000;
B_range = 50;
beta=1;
MIN_PAUSE = 30;
MAX_PAUSE = 3600;
for seed=1:1
    for n_wp=[100]
        for v_Hurst=0.75
%% Generate user movements
            trace=SLAW_MATLAB(dist_alpha, num_user, size_max, n_wp, v_Hurst, ...
                simulationTime / 60, B_range, beta, MIN_PAUSE, MAX_PAUSE);

            filename = strcat('./Heterogeneity/C',num2str(n_wp),'nUser',num2str(num_user));
            save(filename)
        end
    end
end