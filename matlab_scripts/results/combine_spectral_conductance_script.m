%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Livia M. Correa
% Combine spectral conductance 
% Updated 07/11/23 by Livia Correa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
show_axes = 1;
results_path_1 = 'split_1';
results_path_2 = 'split_2';
results_path_3 = '23_07_11_10_51_23';
results_path_4 = '23_07_11_10_51_33';
results_path_5 = '23_07_11_10_52_12';

% Import data
cond_table_1 = readmatrix([results_path_1 '/G_omega_bulk_12_3.000000e+02K.csv']);  
cond_table_2 = readmatrix([results_path_2 '/G_omega_bulk_12_3.000000e+02K.csv']);  
cond_table_3 = readmatrix([results_path_3 '/G_omega_bulk_12_3.000000e+02K.csv']);  
cond_table_4 = readmatrix([results_path_4 '/G_omega_bulk_12_3.000000e+02K.csv']);  
cond_table_5 = readmatrix([results_path_5 '/G_omega_bulk_12_3.000000e+02K.csv']);  
%modes = readtable(['split_modes.xlsx']);

%frequency= table2array(modes(:,1));
frequency_1 = cond_table_1(:,1);
cond_1 = cond_table_1(:,2);

frequency_2 = cond_table_2(:,1);
cond_2 = cond_table_2(:,2);

frequency_3 = cond_table_3(:,1);
cond_3 = cond_table_3(:,2);

frequency_4 = cond_table_4(:,1);
cond_4 = cond_table_4(:,2);

frequency_5 = cond_table_5(:,1);
cond_5 = cond_table_5(:,2);

frequency = [frequency_1;frequency_2;frequency_3;frequency_4;frequency_5];
cond = [cond_1;cond_2;cond_3;cond_4;cond_5];

split_conductance_fig = figure(1);
semilogy(frequency,cond,'-o','LineWidth',1)
xlabel('frequency', 'fontsize', 12)
ylabel('conductance', 'fontsize', 12)


fig_path = ['fig_spectral_conductance.png'];
saveas(split_conductance_fig, fig_path)