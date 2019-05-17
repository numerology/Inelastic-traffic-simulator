% greet_script simulate gain of GREET practical approach under different scenarios
clc, close all, clear all
parpool('local', 40);
warning('off','all');

meanMinFactor = 1:10;
simulationTime = 1000;

poutageGain = zeros(6, length(meanMinFactor));
utilGain = zeros(6, length(meanMinFactor));

for netProfile = 1:6
    fprintf('beginning configuration no. %d\n', netProfile);
    [poutageGain(netProfile, :), utilGain(netProfile, :)] = ...
        GREETsimulation(meanMinFactor, simulationTime, netProfile, 2);
end

%% Plot
datestring = datestr(now, 30);
figure(1);
hold on; grid on;
h1 = plot(meanMinFactor, poutageGain(1, :), 'bo-');
h2 = plot(meanMinFactor, poutageGain(2, :), 'r+-');
h3 = plot(meanMinFactor, poutageGain(3, :), 'kd-');
h4 = plot(meanMinFactor, poutageGain(4, :), 'gh-');
h5 = plot(meanMinFactor, poutageGain(5, :), 'c*-.');
h6 = plot(meanMinFactor, poutageGain(5, :), 'mp:');
set([h1 h2 h3 h4 h5 h6], 'LineWidth', 2, 'MarkerSize', 10);
legend('Scenario 1', 'Scenario 2', 'Scenario 3', 'Scenario 4', 'Scenario 5', ...
    'Scenario 6');
xlabel('mean rate requirement / minimum rate requirement');
ylabel('P(outage) gain over SCPF');
savefig(sprintf('figs/GREET-poutage-multicases-%s.fig', datestring))

figure(2);
hold on; grid on;
h1 = plot(meanMinFactor, utilGain(1, :), 'bo-');
h2 = plot(meanMinFactor, utilGain(2, :), 'r+-');
h3 = plot(meanMinFactor, utilGain(3, :), 'kd-');
h4 = plot(meanMinFactor, utilGain(4, :), 'gh-');
h5 = plot(meanMinFactor, utilGain(5, :), 'c*-.');
h6 = plot(meanMinFactor, utilGain(5, :), 'mp:');
set([h1 h2 h3 h4 h5 h6], 'LineWidth', 2, 'MarkerSize', 10);
legend('Scenario 1', 'Scenario 2', 'Scenario 3', 'Scenario 4', 'Scenario 5', ...
    'Scenario 6');
xlabel('mean rate requirement / minimum rate requirement');
ylabel('utility gain over GPS');
savefig(sprintf('figs/GREET-util-multicases-%s.fig', datestring))
