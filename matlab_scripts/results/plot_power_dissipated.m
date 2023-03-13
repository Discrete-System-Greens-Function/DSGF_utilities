%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Power dissipated using DSGF-c
% Developed by RETL group at the University of Utah, USA
% Last updated 02/03/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
close all;
show_axes = 1;

% mesh=0 for uniform, mesh=1 for non-uniform
mesh = 1;

%Insert the path shown in the terminal.
results_folder = 'results/thin-films_352_subvolumes/SiO2_d_5.00e-07/23_02_24_10_08_42';%'results/thin-films_2_subvolumes/SiO2_d_5.00e-07/23_02_23_11_51_10';%%'results/thin-films_500_subvolumes/SiC_d_1.00e-07/23_02_20_09_54_11';%'results/thin-films_10000_subvolumes/SiC_d_1.00e-07/23_02_17_16_50_46';%'results/thin-films_6000_subvolumes/SiC_d_1.00e-07/23_02_17_12_19_53';%'results/thin-films_4000_subvolumes/SiC_d_1.00e-07/23_02_17_10_49_28';%'results/thin-films_3000_subvolumes/SiC_d_1.00e-07/23_02_17_09_51_50';%'results/thin-films_2000_subvolumes/SiC_d_1.00e-07/23_02_17_09_32_54';%'results/thin-films_1000_subvolumes/SiC_d_1.00e-07/23_02_17_09_13_11';

% Complete directory to import the results
back = cd;
results_path = fullfile(back, '../../', results_folder);

% Import data
R = readmatrix([results_path '/vector_subvolumes_lattice.csv']);  % Make discretized lattice for shape #1 
Q_total_subvol = readmatrix([results_path '/power_dissipated.csv']);
delta_V_vector = readmatrix([results_path '/vector_subvolumes_volume.csv']);
L_sub = delta_V_vector.^(1/3);  

% Set heatmap color axis limits
% To normalize the power dissipated, we can change this limit.
abs_limit = max(abs(Q_total_subvol));

%c_limits = [-abs_limit , abs_limit];
c_limits = [-abs_limit/abs_limit , abs_limit/abs_limit];
%c_limits = [min(Q_total_subvol), max(Q_total_subvol)];

%close all

N_total = size(delta_V_vector,1);
if mesh == 1
i_limit = 1;
for i = 2:N_total
    if(delta_V_vector(i) ~= delta_V_vector(i-1))
        N_subvol_lim(i_limit) = i-1
        i_limit=i_limit+1;
    end
end    

Q_1_c = -Q_total_subvol(1:N_subvol_lim(1))/abs_limit;%/abs_limit
Q_1_r = -Q_total_subvol(N_subvol_lim(1)+1:N_subvol_lim(2))/abs_limit;
Q_2_r = -Q_total_subvol(N_subvol_lim(2)-N_subvol_lim(1)+1:N_subvol_lim(2))/abs_limit;
Q_2_c = -Q_total_subvol(N_subvol_lim(2)+1:N_total)/abs_limit;

R_1_c = R(1:N_subvol_lim(1),:);
R_1_r = R(N_subvol_lim(1)+1:N_subvol_lim(2),:);
R_2_r = R(N_subvol_lim(2)-N_subvol_lim(1)+1:N_subvol_lim(2),:);
R_2_c = R(N_subvol_lim(2)+1:N_total,:);
end

% Subvolume heat map for full particles (VIEW 1)
power_dissipated = figure(1);
if mesh == 0
[vert, fac] = voxel_image( R, L_sub(1), [], [], [], [], 'heatmap', Q_total_subvol/abs_limit, c_limits );
end
if mesh == 1
[vert, fac] = voxel_image( R_1_c, L_sub(1), [], [], [], [], 'heatmap', Q_1_c, c_limits );
[vert, fac] = voxel_image( R_1_r, L_sub(N_subvol_lim(1)+1), [], [], [], [], 'heatmap', Q_1_r, c_limits );
[vert, fac] = voxel_image( R_2_r, L_sub(N_subvol_lim(2)-N_subvol_lim(1)+1), [], [], [], [], 'heatmap',Q_2_r , c_limits );
[vert, fac] = voxel_image( R_2_c, L_sub(N_total), [], [], [], [], 'heatmap', Q_2_c, c_limits );
end
xlabel('x-axis (m)');
ylabel('y-axis (m)');
zlabel('z-axis (m)');
if show_axes == 0
    grid off
    axis off
    colorbar off
end
dim = [.2 .5 .3 .3];
str = ['Q_{Max} = '  num2str(abs_limit) 'W'];
annotation('textbox',dim,'String',str,'FitBoxToText','on');
%view(2)
view(-30,35)

fig_path_power_dissipated = [results_path '/fig_normalized_power_dissipated.fig'];
saveas(power_dissipated, fig_path_power_dissipated)