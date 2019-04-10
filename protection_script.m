clear;
share = 0.5; % share of slice 1.
advIntensityVec = 1:0.2:3;
repetition = 20;
rateVec1Equal = zeros(repetition, length(advIntensityVec));
rateVec2Equal = zeros(repetition, length(advIntensityVec));
rateVec1Dps = zeros(repetition, length(advIntensityVec));
rateVec2Dps = zeros(repetition, length(advIntensityVec));

gcp;

for i = 1:length(advIntensityVec)
    for r = 1:repetition
        in(r) = Simulink.SimulationInput('SCPF_simple');
        in(r) = in(r).setVariable('share1', share);
        in(r) = in(r).setVariable('seed', r);
        in(r) = in(r).setVariable('advArrivalRate', advIntensityVec(i));
        in(r) = in(r).setVariable('weightScheme', 1);
        in(r) = setModelParameter(in(r), 'SaveOutput','on', ...
            'OutputSaveName','out', 'StiffnessThreshold', 20000, 'StopTime', '50');
    end 
    for r = (repetition + 1):(2 * repetition)
        in(r) = Simulink.SimulationInput('SCPF_simple');
        in(r) = in(r).setVariable('share1', share);
        in(r) = in(r).setVariable('seed', r - repetition);
        in(r) = in(r).setVariable('advArrivalRate', advIntensityVec(i));
        in(r) = in(r).setVariable('weightScheme', 2);
        in(r) = setModelParameter(in(r), 'SaveOutput','on', ...
            'OutputSaveName','out', 'StiffnessThreshold', 20000, 'StopTime', '50');
    end
    simOut = parsim(in, 'ShowProgress', 'on');
    
    for r = 1: repetition
        rateVec1Equal(r, i) = simOut(r).out{1}.Values.Data(1, 1, end);
        rateVec2Equal(r, i) = simOut(r).out{1}.Values.Data(1, 2, end);
    end
    for r = (repetition + 1):(2 * repetition)
        rateVec1Dps(r - repetition, i) = simOut(r).out{1}.Values.Data(1, 1, end);
        rateVec2Dps(r - repetition, i) = simOut(r).out{1}.Values.Data(1, 2, end);
    end
    
end

rate1Equal = mean(rateVec1Equal, 1);
rate2Equal = mean(rateVec2Equal, 1);
rate1Dps = mean(rateVec1Dps, 1);
rate2Dps = mean(rateVec2Dps, 1);

datestring = datestr(now, 30);

figure()
hold on
plot(advIntensityVec, rate1Equal, 'b+-');
plot(advIntensityVec, rate2Equal, 'bx-');
plot(advIntensityVec, rate1Dps, 'r+-');
plot(advIntensityVec, rate2Dps, 'rx-');
xlabel('Slice 2 traffic intensity');
ylabel('mean service rate');
legend('SCPF - slice 1', 'SCPF - slice 2', 'DPS - slice1', 'DPS - slice2');
savefig(sprintf('protection-traffic-%s.fig', datestring));
