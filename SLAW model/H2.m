
%% H2
% Fixed
num_user=180*2;
size_max=1000;
Thours=0.002;
MIN_PAUSE=3;
MAX_PAUSE=15;
beta=1;
for seed=17:30
        i
        %% i dependent
        n_wp=50%50+450*i;
        B_range=100%25+250*i;
        dist_alpha=3%1.5+4.5*i;
        v_Hurst=0.55%+0.4*i;
        try
            %% Generate user movements
            trace=SLAW_MATLAB(dist_alpha, num_user, size_max, n_wp, v_Hurst,...
                Thours, B_range, beta, MIN_PAUSE, MAX_PAUSE);
            filename = strcat('./Hetereogeneity/H2_','seed',num2str(seed));
            
            save(filename)
        catch
           disp('bad') 
        end
    
end
