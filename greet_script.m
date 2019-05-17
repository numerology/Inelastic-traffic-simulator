% greet_script simulate gain of GREET practical approach under different scenarios
clc, close all, clear all
parpool('local', 40);
warning('off','all');

meanMinFactor = 1:10;
simulationTime = 1000;

poutageGain = zeros(6, length(meanMinFactor));
poutageError = zeros(6, length(meanMinFactor));
utilGain = zeros(6, length(meanMinFactor));
utilError = zeros(6, length(meanMinFactor));
elasticShares = zeros(6, length(meanMinFactor));

for netProfile = 1:6
    fprintf('beginning configuration no. %d\n', netProfile);
    [poutageGain(netProfile, :), utilGain(netProfile, :), ...
        poutageError(netProfile, :), utilError(netProfile, :), elasticShares(netProfile, :)] = ...
        GREETsimulation(meanMinFactor, simulationTime, netProfile, 2);
end

%% Plot
datestring = datestr(now, 30);
figure(1);
hold on; grid on;
h1 = errorbar(elasticShares(1, :), poutageGain(1, :), poutageError(1, :), 'bo-');
h2 = errorbar(elasticShares(2, :), poutageGain(2, :), poutageError(2, :), 'r+-');
h3 = errorbar(elasticShares(3, :), poutageGain(3, :), poutageError(3, :), 'kd-');
h4 = errorbar(elasticShares(4, :), poutageGain(4, :), poutageError(4, :), 'gh-');
h5 = errorbar(elasticShares(5, :), poutageGain(5, :), poutageError(5, :), 'c*-.');
h6 = errorbar(elasticShares(6, :), poutageGain(6, :), poutageError(6, :), 'mp:');
set([h1 h2 h3 h4 h5 h6], 'LineWidth', 2, 'MarkerSize', 10);
legend('Scenario 1', 'Scenario 2', 'Scenario 3', 'Scenario 4', 'Scenario 5', ...
    'Scenario 6');
xlabel('Overall share of elastic slices');
ylabel('P(outage) gain over share-based approach');
savefig(sprintf('figs/GREET-poutage-multicases-%s.fig', datestring))

figure(2);
hold on; grid on;
h1 = errorbar(elasticShares(1, :), utilGain(1, :), utilError(1, :), 'bo-');
h2 = errorbar(elasticShares(2, :), utilGain(2, :), utilError(2, :), 'r+-');
h3 = errorbar(elasticShares(3, :), utilGain(3, :), utilError(3, :), 'kd-');
h4 = errorbar(elasticShares(4, :), utilGain(4, :), utilError(4, :), 'gh-');
h5 = errorbar(elasticShares(5, :), utilGain(5, :), utilError(5, :), 'c*-.');
h6 = errorbar(elasticShares(6, :), utilGain(6, :), utilError(6, :), 'mp:');
set([h1 h2 h3 h4 h5 h6], 'LineWidth', 2, 'MarkerSize', 10);
legend('Scenario 1', 'Scenario 2', 'Scenario 3', 'Scenario 4', 'Scenario 5', ...
    'Scenario 6');
xlabel('Overall share of elastic slices');
ylabel('utility gain over reservation-based approach');
savefig(sprintf('figs/GREET-util-multicases-%s.fig', datestring))
