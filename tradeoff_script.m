% script to get the tradeoff comparison among the 4 benchmarks.
clc, close all, clear all
parpool('local', 40);
warning('off','all');

satVec = [1 3 5 7 10 13 15]; % subject to what type of SLAW model to use
simulationTime = 3000;

[poutage, utility] = GREETsimulation(satVec, simulationTime, 2, 1);

%% Plot
datestring = datestr(now, 30);
figure()
hold on; grid on;
h1 = plot(utility(1, :), 1 - poutage(1, :), 'ro-');
h2 = plot(utility(2, :), 1 - poutage(2, :), 'ch-');
h3 = plot(utility(3, :), 1 - poutage(3, :), 'b+-');
h4 = plot(utility(4, :), 1 - poutage(4, :), 'gd-');
legend('GREET-PRA', 'GREET-OPT', 'SCPF', 'GPS');
set([h1 h2 h3 h4], 'LineWidth', 2, 'MarkerSize', 10);
xlabel('Utility');
ylabel('1 - P(outage)');
savefig(sprintf('figs/GREET-poutage-utility-tradeoff-%s.fig', datestring))