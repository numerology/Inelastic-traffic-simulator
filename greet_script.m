% greet_script simulate gain of GREET practical approach under different scenarios
clc, close all, clear all
%parpool('local', 40);
warning('off','all');

meanMinFactor = 1:10;
% Simulation time span, cannot exceed 1000.
simulationTime = 1000;

poutageGainAgainstSCPF = zeros(6, length(meanMinFactor));
poutageErrorAgainstSCPF = zeros(6, length(meanMinFactor));
poutageGainAgainstReservation = zeros(6, length(meanMinFactor));
poutageErrorAgainstReservation = zeros(6, length(meanMinFactor));
utilGainAgainstReservation = zeros(6, length(meanMinFactor));
utilErrorAgainstReservation = zeros(6, length(meanMinFactor));
utilGainAgainstSCPF = zeros(6, length(meanMinFactor));
utilErrorAgainstSCPF = zeros(6, length(meanMinFactor));
elasticShares = zeros(6, length(meanMinFactor));

% Network profile 1, 3, 4, 6 are what we used in our WiOpt'20 submission.
% Please refer to GREETsimulation re: what they really mean.
for netProfile = [1 3 4 6]
    fprintf('beginning configuration no. %d\n', netProfile);
    [poutageGainAgainstSCPF(netProfile, :), ...
        poutageGainAgainstReservation(netProfile, :), ...
        utilGainAgainstReservation(netProfile, :), ...
        utilGainAgainstSCPF(netProfile, :), ...
        poutageErrorAgainstSCPF(netProfile, :), ...
        poutageErrorAgainstReservation(netProfile, :), ...
        utilErrorAgainstReservation(netProfile, :), ...
        utilErrorAgainstSCPF(netProfile, :), ...
        elasticShares(netProfile, :)] = ...
        GREETsimulation(meanMinFactor, simulationTime, netProfile, 2);
end

%% New plot
datestring = datestr(now, 30);
figure(3);
hold on; grid on;
h1 = errorbar(elasticShares(1, :), poutageGainAgainstReservation(1, :), poutageErrorAgainstReservation(1, :), 'bo-');
% h2 = errorbar(elasticShares(2, :), poutageGainAgainstSCPF(2, :), poutageErrorAgainstSCPF(2, :), 'r+-');
h3 = errorbar(elasticShares(3, :), poutageGainAgainstReservation(3, :), poutageErrorAgainstReservation(3, :), 'kd-');
h4 = errorbar(elasticShares(4, :), poutageGainAgainstReservation(4, :), poutageErrorAgainstReservation(4, :), 'gh-');
% h5 = errorbar(elasticShares(5, :), poutageGainAgainstSCPF(5, :), poutageErrorAgainstSCPF(5, :), 'c*-.');
h6 = errorbar(elasticShares(6, :), poutageGainAgainstReservation(6, :), poutageErrorAgainstReservation(6, :), 'mp:');
set([h1 h3 h4 h6], 'LineWidth', 2, 'MarkerSize', 10);
legend('Uniform', 'Het Aligned', 'Het Orthogonal', ...
    'Mixed');
xlabel('Overall share of elastic slices');
ylabel('P(outage) gain over reservation-based approach');
xlim([0 30])
title('Gain in P(outage) over reservation-based approach')
savefig(sprintf('figs/GREET-poutage-multicases-new-%s.fig', datestring))

figure(4);
hold on; grid on;
h1 = errorbar(elasticShares(1, :), utilGainAgainstSCPF(1, :), utilErrorAgainstSCPF(1, :), 'bo-');
% h2 = errorbar(elasticShares(2, :), utilGainAgainstReservation(2, :), utilErrorAgainstReservation(2, :), 'r+-');
h3 = errorbar(elasticShares(3, :), utilGainAgainstSCPF(3, :), utilErrorAgainstSCPF(3, :), 'kd-');
h4 = errorbar(elasticShares(4, :), utilGainAgainstSCPF(4, :), utilErrorAgainstSCPF(4, :), 'gh-');
% h5 = errorbar(elasticShares(5, :), utilGainAgainstReservation(5, :), utilErrorAgainstReservation(5, :), 'c*-.');
h6 = errorbar(elasticShares(6, :), utilGainAgainstSCPF(6, :), utilErrorAgainstSCPF(6, :), 'mp:');
set([h1 h3 h4 h6], 'LineWidth', 2, 'MarkerSize', 10);
legend('Uniform', 'Het Aligned', 'Het Orthogonal', ...
    'Mixed');
xlabel('Overall share of elastic slices');
ylabel('utility gain over share-based approach');
xlim([0 30])
title('Gain in utility over share-based approach')
savefig(sprintf('figs/GREET-util-multicases-new-%s.fig', datestring))

%% Plot
datestring = datestr(now, 30);
figure(1);
hold on; grid on;
h1 = errorbar(elasticShares(1, :), poutageGainAgainstSCPF(1, :), poutageErrorAgainstSCPF(1, :), 'bo-');
h2 = errorbar(elasticShares(2, :), poutageGainAgainstSCPF(2, :), poutageErrorAgainstSCPF(2, :), 'r+-');
h3 = errorbar(elasticShares(3, :), poutageGainAgainstSCPF(3, :), poutageErrorAgainstSCPF(3, :), 'kd-');
h4 = errorbar(elasticShares(4, :), poutageGainAgainstSCPF(4, :), poutageErrorAgainstSCPF(4, :), 'gh-');
h5 = errorbar(elasticShares(5, :), poutageGainAgainstSCPF(5, :), poutageErrorAgainstSCPF(5, :), 'c*-.');
h6 = errorbar(elasticShares(6, :), poutageGainAgainstSCPF(6, :), poutageErrorAgainstSCPF(6, :), 'mp:');
set([h1 h2 h3 h4 h5 h6], 'LineWidth', 2, 'MarkerSize', 10);
legend('Scenario 1', 'Scenario 2', 'Scenario 3', 'Scenario 4', 'Scenario 5', ...
    'Scenario 6');
xlabel('Overall share of elastic slices');
ylabel('P(outage) gain over share-based approach');
savefig(sprintf('figs/GREET-poutage-multicases-%s.fig', datestring))

figure(2);
hold on; grid on;
h1 = errorbar(elasticShares(1, :), utilGainAgainstReservation(1, :), utilErrorAgainstReservation(1, :), 'bo-');
h2 = errorbar(elasticShares(2, :), utilGainAgainstReservation(2, :), utilErrorAgainstReservation(2, :), 'r+-');
h3 = errorbar(elasticShares(3, :), utilGainAgainstReservation(3, :), utilErrorAgainstReservation(3, :), 'kd-');
h4 = errorbar(elasticShares(4, :), utilGainAgainstReservation(4, :), utilErrorAgainstReservation(4, :), 'gh-');
h5 = errorbar(elasticShares(5, :), utilGainAgainstReservation(5, :), utilErrorAgainstReservation(5, :), 'c*-.');
h6 = errorbar(elasticShares(6, :), utilGainAgainstReservation(6, :), utilErrorAgainstReservation(6, :), 'mp:');
set([h1 h2 h3 h4 h5 h6], 'LineWidth', 2, 'MarkerSize', 10);
legend('Scenario 1', 'Scenario 2', 'Scenario 3', 'Scenario 4', 'Scenario 5', ...
    'Scenario 6');
xlabel('Overall share of elastic slices');
ylabel('utility gain over reservation-based approach');
savefig(sprintf('figs/GREET-util-multicases-%s.fig', datestring))
