%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spectral conductance using DSGF-c
% Developed by RETL group at the University of Utah, USA
% Last updated 01/18/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
show_axes = 1;
% Directory where results files will be saved (string)
back = cd;

%Insert the path shown in the terminal.
results_folder = 'results/thin-films_352_subvolumes/SiO2_d_5.00e-07/23_02_24_10_08_42';%'results/thin-films_10000_subvolumes/SiC_d_1.00e-07/23_02_17_16_50_46';%'results/thin-films_6000_subvolumes/SiC_d_1.00e-07/23_02_17_12_19_53';%'results/thin-films_4000_subvolumes/SiC_d_1.00e-07/23_02_17_10_49_28';%'results/thin-films_3000_subvolumes/SiC_d_1.00e-07/23_02_17_09_51_50';%'results/thin-films_2000_subvolumes/SiC_d_1.00e-07/23_02_17_09_32_54';%'results/thin-films_1000_subvolumes/SiC_d_1.00e-07/23_02_17_09_23_17';%'results/thin-films_1000_subvolumes/SiC_d_1.00e-07/23_02_17_09_13_11';

% Complete directory to import the results
results_path = fullfile(back, '../../', results_folder);

% Import data
temperature = readmatrix([results_path '/total_conductance.csv'],'Range','A:A');  % results_path '/spectral_conductance_3.000000e+02K.csv'
G_12 = readmatrix([results_path '/total_conductance.csv'],'Range','B:B');

% Plot spectral SGF results
G_total_12 = figure(1);
loglog(temperature, G_12,'-o', 'linewidth', 2, 'markersize', 5)%'-o'
xlabel('Temperature, [K]')
ylabel('Conductance [W/K]')
xticks([200,250,300,350,400])
%axis([200 400 0.4e-11 3e-11])
set(gca, 'fontsize', 16)
grid on
hold off

fig_path_total_conductance = [results_path '/fig_total_conductance.fig']; %
saveas(G_total_12, fig_path_total_conductance) %