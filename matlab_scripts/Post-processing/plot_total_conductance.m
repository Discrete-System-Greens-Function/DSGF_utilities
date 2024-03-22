%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spectral conductance using DSGF-c
% Developed by RETL group at the University of Utah, USA
% Last updated 03/22/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
show_axes = 1;

%Insert the path shown in the terminal.
results_path = 'results/user_defined_576_subvolumes/Si3N4_d_5.00e-07/23_11_20_08_51_47';

%              'results/sample_2_subvolumes/SiC_d_1.00e-07/23_11_15_10_12_24';
%              'results/sample_416_subvolumes/SiO2_d_1.50e-07/23_11_16_08_38_31';
%              'results/sample_2000_subvolumes/SiO2_d_1.00e-07/23_11_16_08_59_44';
%              'results/user_defined_576_subvolumes/Si3N4_d_5.00e-07/23_11_20_08_51_47';

% Complete directory to import the results
% Directory where results files will be saved (string)
%back = cd;
%results_path = fullfile(back, '../../', results_folder);

% Import data
temperature = readmatrix([results_path '/G_t_AB.csv'],'Range','A:A');  % results_path '/spectral_conductance_3.000000e+02K.csv'
G_12 = readmatrix([results_path '/G_t_AB.csv'],'Range','B:B');

% Plot spectral SGF results
G_total_12 = figure(1);
loglog(temperature, G_12,'-o', 'linewidth', 2, 'markersize', 5)%'-o'
xlabel('T [K]')
ylabel('G_{t,AB} [W/K]')
xticks([200,250,300,350,400])
%axis([200 400 0.4e-11 3e-11])
set(gca, 'fontsize', 16)
grid on
hold off

fig_path_total_conductance = [results_path '/fig_total_conductance.fig']; %
saveas(G_total_12, fig_path_total_conductance) %