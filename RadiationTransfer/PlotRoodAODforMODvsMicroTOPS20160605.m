%% Plot MODTRAN AOTs on 2016-06-5 vs MicroTOPS
clear all
close all
% Plotting sizes
theFontSize = 14;
theLineWidth = 2.0;
% Data saved with
% save(['RoodeplaatAOT' MeasureDate IHAZEModel '.mat'], 'AODWv', 'Diff2GlobAOD', 'AOTwv', 'AOT');
load RoodeplaatAOT20160605Tuned.mat
Diff2GlobAODTuned = Diff2GlobAOD;
load RoodeplaatAOT20160605Urban.mat
Diff2GlobAODUrban = Diff2GlobAOD;
load RoodeplaatAOT20160605Rural.mat
Diff2GlobAODRural = Diff2GlobAOD;


figure;
plot(AODWv, Diff2GlobAODTuned, '--', AODWv, Diff2GlobAODRural, AODWv, Diff2GlobAODUrban, ':', AOTwv, AOT, 'or', 'LineWidth', theLineWidth);
set(gca, 'LineWidth', theLineWidth);
set(gca, 'FontSize', theFontSize);
title(['Vertical AOD, Roodeplaat 2016-06-05'])
xlabel('Wavelength [nm]');
ylabel('AOD')
legend('MODTRAN Tuned', 'MODTRAN Rural', 'MODTRAN Urban', 'MicroTOPS', 'location', 'best');
axis([350, 900, 0, 1])
grid();