thresholdVec = 0.1:0.1:1;
slice1Perf = zeros(size(thresholdVec));
slice2Perf = zeros(size(thresholdVec));
slice1Gain = zeros(size(thresholdVec));
slice2Gain = zeros(size(thresholdVec));
for i = 1:size(thresholdVec, 2)
    fprintf('current threshold %f\n', thresholdVec(i));
    parfor t = 1:simulationTime
        [r,f,b] = GPS(NetSettings, OpSettings, [capacityPerUser(:,t)]', ...
            [bs(:,t)]');
        rates_GPS(:,t)=r;
        fractions_GPS(:,t)=f;
        btd_GPS(:,t)=b;
        [r,f,b] = newsharing(NetSettings, OpSettings, [capacityPerUser(:,t)]',...
        [bs(:,t)]', thresholdVec(i));
        rates_new(:,t) = r;
        fractions_new(:,t) = f;
        btd_new(:,t) = b;
    end
    
    slice1Perf(i) = mean(mean(btd_new(OpSettings.ops_belongs == 1, :, :)));
    slice2Perf(i) = mean(mean(btd_new(OpSettings.ops_belongs == 2, :, :)));
    slice1Gain(i) = mean(mean(btd_GPS(OpSettings.ops_belongs == 1, :, :)))...
        / slice1Perf(i);
    slice2Gain(i) = mean(mean(btd_GPS(OpSettings.ops_belongs == 2, :, :)))...
        / slice2Perf(i);
end

figure()
hold on
plot(thresholdVec, slice1Perf);
plot(thresholdVec, slice2Perf);
title('BTD of new sharing')
legend('slice1 - uniform', 'slice2 - het');

figure()
hold on
plot(thresholdVec, slice1Gain);
plot(thresholdVec, slice2Gain);
title('Gain vs. GPS of new sharing')
legend('slice1 - uniform', 'slice2 - het');