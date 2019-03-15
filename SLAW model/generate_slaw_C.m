% Generate clustered user movement trace
%% Users settings
dist_alpha=6;
% need to support at most 19 * 30
num_user = 190 * 3;
size_max = 1000;
simulationTime = 5000;
B_range = 300;
beta=1;
MIN_PAUSE = 100;
MAX_PAUSE = 300;
for seed=1:1
    for n_wp=[100]
        for v_Hurst=0.9
%% Generate user movements
            trace=SLAW_MATLAB(dist_alpha, num_user, size_max, n_wp, v_Hurst, ...
                simulationTime / 60, B_range, beta, MIN_PAUSE, MAX_PAUSE);

            filename = strcat('./Heterogeneity/v90alpha',num2str(dist_alpha),'nUser',num2str(num_user));
            save(filename)
        end
    end
end