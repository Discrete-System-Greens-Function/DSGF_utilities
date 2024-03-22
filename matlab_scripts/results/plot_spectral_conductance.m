%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spectral conductance using DSGF-c
% Developed by RETL group at the University of Utah, USA
% Last updated 01/18/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
close all;
show_axes = 1;

%Insert the path shown in the terminal.
results_folder = 'results/thin-films_352_subvolumes/SiO2_d_5.00e-07/23_02_24_10_08_42';

% Complete directory to import the results
back = cd;
results_path = fullfile(back, '../../', results_folder);

% Import data
frequency = readmatrix([results_path '/spectral_conductance_3.000000e+02K.csv'],'Range','A:A');  % results_path '/spectral_conductance_3.000000e+02K.csv'
G_12_omega_200 = readmatrix([results_path '/spectral_conductance_2.000000e+02K.csv'],'Range','B:B');
G_12_omega_250 = readmatrix([results_path '/spectral_conductance_2.500000e+02K.csv'],'Range','B:B');
G_12_omega_300 = readmatrix([results_path '/spectral_conductance_3.000000e+02K.csv'],'Range','B:B');
G_12_omega_350 = readmatrix([results_path '/spectral_conductance_3.500000e+02K.csv'],'Range','B:B');
G_12_omega_400 = readmatrix([results_path '/spectral_conductance_4.000000e+02K.csv'],'Range','B:B');

% Plot spectral SGF results
G_spectral_12 = figure(1);
loglog(frequency, G_12_omega_200,'-', 'linewidth', 2, 'markersize', 3)%'-o'
%loglog(omega_SGF*(10^6), G_12_omega_SGF.*(10^9).*(1e6), '-', 'linewidth', 2, 'markersize', 5)
%xlabel('Wavelength, \lambda [um]')
%ylabel('Spectral conductance, G_1_-_>_2(\lambda) [nWK^-^1um^-^1]')
hold on
loglog(frequency, G_12_omega_250,'-', 'linewidth', 2, 'markersize', 3) % use 'o' for parallel and '-' for serial
loglog(frequency, G_12_omega_300,'-', 'linewidth', 2, 'markersize', 3)
loglog(frequency, G_12_omega_350,'-', 'linewidth', 2, 'markersize', 3)
loglog(frequency, G_12_omega_400,'-', 'linewidth', 2, 'markersize', 3)
xlabel('Frequency, \omega [rad/s]')
ylabel('G_1_-_>_2(\omega) [W/K]')
title('Spectral conductance for different temperatures ')
legend('200 K','250 K','300 K','350 K','400 K', 'location', 'best');
set(gca, 'fontsize', 16)
grid on
hold off

fig_path_spectral_conductance = [results_path '/fig_spectral_conductance.fig']; %
saveas(G_spectral_12, fig_path_spectral_conductance) %
