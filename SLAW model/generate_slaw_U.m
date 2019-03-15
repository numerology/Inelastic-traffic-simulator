%% Users settings
dist_alpha=0.5;
num_user= 190 * 3;
size_max=1000;
Thours = 5000 / 60;
B_range=1;
beta=1;
MIN_PAUSE=1;
MAX_PAUSE=6;
for seed=1:1
    for n_wp=[100]%75 100 125 150 200 400]
        for v_Hurst=0.53;
%% Generate user movements
trace=SLAW_MATLAB(dist_alpha, num_user, size_max, n_wp, v_Hurst,Thours, B_range, beta, MIN_PAUSE, MAX_PAUSE);

filename = strcat('./Heterogeneity/U',num2str(n_wp),'seed',num2str(seed));
save(filename)
        end
    end
end