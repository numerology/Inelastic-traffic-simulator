

%% Fixed
num_user=180*2;
size_max=1000;
Thours=0.002;
MIN_PAUSE=3;
MAX_PAUSE=15;
beta=1;

for seed=1%8:10
    counter=0;
    for i=1%0.1:0.1:0.7
        i
        counter=counter+1;
        %% i dependent
        n_wp=75%50+450*i;
        B_range=60%25+250*i;
        dist_alpha=3.5%1.5+4.5*i;
        v_Hurst=0.95%0.55+0.4*i;
        subplot(3,3,counter)
        try
        %% Generate user movements
        trace=SLAW_MATLAB(dist_alpha, num_user, size_max, n_wp, v_Hurst,...
            Thours, B_range, beta, MIN_PAUSE, MAX_PAUSE);
        filename = strcat('./Hetereogeneity/test',num2str(1))%,'seed',num2str(seed));
        hold on
        plot(trace(:,300,1),trace(:,300,2),'r*')
        save(filename)
        catch
           disp('bad') 
        end
    end
end