%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Power dissipated using DSGF-c
% Developed by RETL group at the University of Utah, USA
% Last updated 02/03/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
close all;
show_axes = 0;

% mesh=0 for sample with same subvolume size, mesh=, for sample with different subvolume size per thermal object, mesh=2 for user-defined with non-uniform discretization , 
mesh = 2;

% Norm =0 for full range, norm=1 for normalized;
norm = 0;

%Insert the path shown in the terminal.
results_path = 'results/user_defined_576_subvolumes/Si3N4_d_5.00e-07/23_11_20_08_51_47';

% Examples:

% 'results/sample_2_subvolumes/SiC_d_1.00e-07/23_11_15_10_12_24'; % use mesh = 0             
% 'results/sample_2000_subvolumes/SiO2_d_1.00e-07/23_11_16_08_59_44'; % use mesh = 0 
% 'results/sample_416_subvolumes/SiO2_d_1.50e-07/23_11_16_08_38_31'; % use mesh = 1 
% 'results/user_defined_576_subvolumes/Si3N4_d_5.00e-07/23_11_20_08_51_47'; % use mesh = 2 

% Complete directory to import the results
%back = cd;
%results_path = fullfile(back, results_folder);

% Import data
R = readmatrix([results_path '/vector_subvolumes_lattice.csv']);  % Make discretized lattice for shape #1 
Q_total_subvol = readmatrix([results_path '/Q_density_subvol.csv']);
delta_V_vector = readmatrix([results_path '/vector_subvolumes_volume.csv']);
L_sub = delta_V_vector.^(1/3);  

% Set heatmap color axis limits
% To normalize the power dissipated, we can change this limit.
abs_limit = max(abs(Q_total_subvol));

if norm == 0 max_limit = 1 ; end
if norm == 1 max_limit = abs_limit; end

c_limits = [-abs_limit/max_limit , abs_limit/max_limit];

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

Q_1_c = -Q_total_subvol(1:N_subvol_lim(1))/max_limit;%/abs_limit
Q_1_r = -Q_total_subvol(N_subvol_lim(1)+1:N_total)/max_limit;

R_1_c = R(1:N_subvol_lim(1),:);
R_1_r = R(N_subvol_lim(1)+1:N_total,:);
end

if mesh == 2
i_limit = 1;
for i = 2:N_total
    if(delta_V_vector(i) ~= delta_V_vector(i-1))
        N_subvol_lim(i_limit) = i-1
        i_limit=i_limit+1;
    end
end    

Q_1_c = -Q_total_subvol(1:N_subvol_lim(1))/max_limit;%/abs_limit
Q_1_r = -Q_total_subvol(N_subvol_lim(1)+1:N_subvol_lim(2))/max_limit;
Q_2_r = -Q_total_subvol(N_subvol_lim(2)-N_subvol_lim(1)+1:N_subvol_lim(2))/max_limit;
Q_2_c = -Q_total_subvol(N_subvol_lim(2)+1:N_total)/max_limit;

R_1_c = R(1:N_subvol_lim(1),:);
R_1_r = R(N_subvol_lim(1)+1:N_subvol_lim(2),:);
R_2_r = R(N_subvol_lim(2)-N_subvol_lim(1)+1:N_subvol_lim(2),:);
R_2_c = R(N_subvol_lim(2)+1:N_total,:);
end

% Subvolume heat map for full particles (VIEW 1)
power_dissipated = figure(1);
if mesh == 0
[vert, fac] = voxel_density( R, L_sub(1), [], [], [], [], 'heatmap', Q_total_subvol/max_limit, c_limits );
end
if mesh == 1
[vert, fac] = voxel_power( R_1_c, L_sub(1), [], [], [], [], 'heatmap', Q_1_c/max_limit, c_limits );
[vert, fac] = voxel_power( R_1_r, L_sub(N_subvol_lim(1)+1), [], [], [], [], 'heatmap', Q_1_r/max_limit, c_limits );
end
if mesh == 2
[vert, fac] = voxel_density( R_1_c, L_sub(1), [], [], [], [], 'heatmap', Q_1_c/max_limit, c_limits );
[vert, fac] = voxel_density( R_1_r, L_sub(N_subvol_lim(1)+1), [], [], [], [], 'heatmap', Q_1_r/max_limit, c_limits );
[vert, fac] = voxel_density( R_2_r, L_sub(N_subvol_lim(2)-N_subvol_lim(1)+1), [], [], [], [], 'heatmap',Q_2_r/max_limit , c_limits );
[vert, fac] = voxel_density( R_2_c, L_sub(N_total), [], [], [], [], 'heatmap', Q_2_c/max_limit, c_limits );
end
xlabel('x-axis (m)');
ylabel('y-axis (m)');
zlabel('z-axis (m)');
if show_axes == 0
    grid off
    axis off
    %colorbar off
end
dim = [.2 .5 .3 .3];
str = ['Q_{Max} = '  num2str(abs_limit) 'W'];
annotation('textbox',dim,'String',str,'FitBoxToText','on');
%view(2)
view(-30,35)

fig_path_power_dissipated = [results_path '/fig_power_density.fig'];
saveas(power_dissipated, fig_path_power_dissipated)