%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lindsay Walter and Livia Correa
% Thin film discretization
% Updated 07/03/24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear, clc, close all
format long
profile on
disp(['Running MATLAB script ' mfilename])

% This code creates discretizations for a system of two cuboids
% separated by a gap.  Results are output in either .txt format or in an
% Excel file based on user input.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%******************************USER INPUTS********************************%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set dimensions and mesh refinement %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

description = '2_cubes';  % Description of geometry

d = 500e-9;     % Distance between the two cuboids [m]
Lx = 500e-9;    % Dimension in x-direction for both the cuboids [m]
Ly = 500e-9;    % Dimension in y-direction for both the cuboids [m]
Lz = 500e-9;    % Dimension in z-direction for both the cuboids [m]

% To generate uniform and nonuniform discretization: 
    % Uniform: set refinements = 0 and n_0 with the number of subvolumes across the smallest dimension given.
    % Nonuniform with 1 refinement: set refinements = 1, n_0 with the number of subvolumes for the less refined mesh, Lx_ref_1 with the more refined length closer to the gap, and n_1 with the number of subvolumes for the more refined mesh. 
    % Nonuniform with 2 refinements: set refinements = 2, n_0 with the number of subvolumes for the less refined mesh, Lx_ref_1 with the intermediate refined length, n_1 with the number of subvolumes for the intermediate refined mesh, Lx_ref_2 with the more refined length closer to the gap, and n_2 with the number of subvolumes for the more refined mesh. 

refinements = 1;    % Choose between 0, 1, or 2. 
n_0 = 4;            % Number of subvolumes for no refinement across the smallest dimension given 
Lx_ref_1 = 250e-9;  % Refinement length if refinement = 1 or 2
n_1 = 8;            % Number of subvolumes  for refinement 1 across the smallest dimension given
Lx_ref_2 = 250e-9;  % Refinement length if refinement = 2 
n_2 = 4;            % Number of subvolumes for refinement 2 across the smallest dimension given 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set type of discretization export %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Save discretization in an Excel file? (0 = "no", 1 = "yes")
save_Excel = 0;

% Save discretization in an .txt file? (0 = "no", 1 = "yes")
save_txt = 1;

% Save figures? (0 = "no", 1 = "yes")
save_fig = 1;


%%%%%%%%%%%%%%%%%%
% Figure options %
%%%%%%%%%%%%%%%%%%

% Show figure axes? (0 = "no", 1 = "yes")
show_axes = 1;

% Show figure titles? (0 = "no", 1 = "yes")
show_titles = 1;

%****************************END USER INPUTS******************************%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cord = 'min'; % String of coordinate along which direction mesh is to be refined. 'Min' refers to the smallest dimension.

if refinements == 0
    % Thermal source #1
    Lx_1 = Lx;         % Length of film #1 
    Ly_1 = Ly;         % Width of film #1
    Lz_1 = Lz;         % Thickness of film #1 
    origin_1 = [0,0,0];    % Point of back, left, bottom of cube [x,y,z]
    cord_1 = cord;          % String of coordinate along which direction mesh is to be refined 
    mesh_1 = n_0;            % Number of subvolumes across the dimension given in the input 'cord' 

    % Thermal source #2
    Lx_2 = Lx_1;              % Length of film #2
    Ly_2 = Ly_1;              % Width of film #2
    Lz_2 = Lz_1;              % Thickness of film #2
    origin_2 = [Lx_1+d,0,0];  % Point of back, left, bottom of cube [x,y,z]
    cord_2 = cord;             % String of coordinate along which direction mesh is to be refined 
    mesh_2 = mesh_1;               % Number of subvolumes across the dimension given in the input 'cord'   
end

if refinements == 1
    % Thermal source #1
    Lx_1_c = Lx-Lx_ref_1;         % Length of film #1 
    Lx_1_r = Lx_ref_1;         % Length of film #1 
    Ly_1 = Ly;         % Width of film #1
    Lz_1 = Lz;         % Thickness of film #1 
    origin_1 = [0,0,0];    % Point of back, left, bottom of cube [x,y,z]
    origin_1_r = [Lx_1_c,0,0];    % Point of back, left, bottom of cube [x,y,z]
    cord_1 = cord;          % String of coordinate along which direction mesh is to be refined 
    mesh_1_c = n_0;            % Number of subvolumes across the dimension given in the input 'cord' 
    mesh_1_r = n_1;            % Number of subvolumes across the dimension given in the input 'cord' 

    % Thermal source #2
    Lx_2_c = Lx_1_c;              % Length of film #2
    Lx_2_r = Lx_1_r;              % Length of film #2
    Ly_2 = Ly_1;              % Width of film #2
    Lz_2 = Lz_1;              % Thickness of film #2
    origin_2_r = [d+Lx_1_c+Lx_1_r,0,0];  % Point of back, left, bottom of cube [x,y,z
    origin_2_c = [d+Lx_1_c+Lx_1_r+Lx_2_r,0,0];  % Point of back, left, bottom of cube [x,y,z]
    cord_2 = cord;             % String of coordinate along which direction mesh is to be refined 
    mesh_2_c = mesh_1_c;               % Number of subvolumes across the dimension given in the input 'cord' 
    mesh_2_r = mesh_1_r;               % Number of subvolumes across the dimension given in the input 'cord'
end

if refinements == 2
    % Thermal source #1
    Lx_1_c = Lx-Lx_ref_1-Lx_ref_2;         % Length of film #1 
    Lx_1_r = Lx_ref_1;         % Length of film #1 
    Lx_1_r_2 = Lx_ref_2;         % Length of film #1 
    Ly_1 = Ly;         % Width of film #1
    Lz_1 = Lz;         % Thickness of film #1 
    origin_1 = [0,0,0];    % Point of back, left, bottom of cube [x,y,z]
    origin_1_r = [Lx_1_c,0,0];    % Point of back, left, bottom of cube [x,y,z]
    origin_1_r_2 = [Lx_1_c+Lx_1_r,0,0];    % Point of back, left, bottom of cube [x,y,z]
    cord_1 = cord;          % String of coordinate along which direction mesh is to be refined 
    mesh_1_c = n_0;            % Number of subvolumes across the dimension given in the input 'cord' 
    mesh_1_r = n_1;            % Number of subvolumes across the dimension given in the input 'cord' 
    mesh_1_r_2 = n_2;            % Number of subvolumes across the dimension given in the input 'cord' 

    % Thermal source #2
    Lx_2_c = Lx_1_c;              % Length of film #2
    Lx_2_r = Lx_1_r;              % Length of film #2
    Lx_2_r_2 = Lx_1_r_2;              % Length of film #2
    Ly_2 = Ly_1;              % Width of film #2
    Lz_2 = Lz_1;              % Thickness of film #2
    origin_2_r_2 = [d+Lx_1_c+Lx_1_r+Lx_1_r_2,0,0];    % Point of back, left, bottom of cube [x,y,z]
    origin_2_r = [d+Lx_1_c+Lx_1_r+Lx_1_r_2+Lx_2_r_2,0,0];  % Point of back, left, bottom of cube [x,y,z
    origin_2_c = [d+Lx_1_c+Lx_1_r+Lx_2_r+Lx_1_r_2+Lx_2_r_2,0,0];  % Point of back, left, bottom of cube [x,y,z]
    cord_2 = cord;             % String of coordinate along which direction mesh is to be refined 
    mesh_2_c = mesh_1_c;               % Number of subvolumes across the dimension given in the input 'cord' 
    mesh_2_r = mesh_1_r;               % Number of subvolumes across the dimension given in the input 'cord'
    mesh_2_r_2 = mesh_1_r_2;               % Number of subvolumes across the dimension given in the input 'cord'
end


%%%%%%%%%%%%%%%%%%%%%%%%%
% Create discretization %
%%%%%%%%%%%%%%%%%%%%%%%%%

if refinements == 0  
    % Combine inputs
    L_1 = [Lx_1, Ly_1, Lz_1]; % Vector containing length of each side of film #1 [m]
    L_2 = [Lx_1, Ly_1, Lz_1]; % Vector containing length of each side of film #1 [m]
    
    [ R_1, delta_V_1, N_1 ] = discretize_cuboid(L_1, origin_1, mesh_1, cord_1);  % Create discretization for cuboid #1
    [ R_2, delta_V_2, N_2 ] = discretize_cuboid(L_2, origin_2, mesh_2, cord_2); % Create discretization for cuboid #2
    r = [R_1; R_2]; % Combine discretization for both films in one matrix
    N = N_1 + N_2; % Total number of subvolumes in both films
    ind_bulk = [1, N_1+1]; % Bulk object start index
    delta_V_vector = [delta_V_1, delta_V_2]; % Vector of volume of each subvolume [m^3]
end   

if refinements == 1
    % Combine inputs
    L_1_c = [Lx_1_c, Ly_1, Lz_1]; % Vector containing length of each side of film #1 [m]
    L_1_r = [Lx_1_r, Ly_1, Lz_1]; % Vector containing length of each side of film #1 [m]
    L_2_r = [Lx_2_r, Ly_2, Lz_2]; % Vector containing length of each side of film #1 [m]
    L_2_c = [Lx_2_c, Ly_2, Lz_2]; % Vector containing length of each side of film #1 [m]
    
    [ R_1_c, delta_V_1_c, N_1_c ] = discretize_cuboid(L_1_c, origin_1, mesh_1_c, cord_1); % Create discretization for coarse part of cuboid #1
    [ R_1_r, R_2_r, delta_V_r, N_r ] = mirror_nonuniform(L_1_r, origin_1_r, mesh_1_r, cord_1, d,Lx,refinements); % Create discretization for mirror the refined parts of the cuboids
    [ R_2_c, delta_V_2_c, N_2_c ] = discretize_cuboid(L_2_c, origin_2_c, mesh_2_c, cord_2); % Create discretization for coarse part of cuboid #2
    r = [R_1_c;R_1_r; R_2_r; R_2_c]; % Combine discretization in one matrix
    N = N_1_c + N_r + N_2_c; % Total number of subvolumes 
    N_1 = N_1_c + N_r/2;
    delta_V_vector = [delta_V_1_c, delta_V_r, delta_V_2_c]; % Vector of volume of each subvolume [m^3]
    
end    

if refinements == 2
    % Combine inputs
    L_1_c = [Lx_1_c, Ly_1, Lz_1]; % Vector containing length of each side of film #1 [m]
    L_1_r = [Lx_1_r, Ly_1, Lz_1]; % Vector containing length of each side of film #1 [m]
    L_1_r_2 = [Lx_1_r_2, Ly_1, Lz_1]; % Vector containing length of each side of film #1 [m]
    L_2_r_2 = [Lx_2_r_2, Ly_2, Lz_2]; % Vector containing length of each side of film #1 [m]
    L_2_r = [Lx_2_r, Ly_2, Lz_2]; % Vector containing length of each side of film #1 [m]
    L_2_c = [Lx_2_c, Ly_2, Lz_2]; % Vector containing length of each side of film #1 [m]
    
    [ R_1_c, delta_V_1_c, N_1_c ] = discretize_cuboid(L_1_c, origin_1, mesh_1_c, cord_1); % Create discretization for coarse part of cuboid #1
    [ R_1_r, delta_V_1_r, N_1_r ] = discretize_cuboid(L_1_r, origin_1_r, mesh_1_r, cord_1); % Create discretization for coarse part of cuboid #1
    [ R_1_r_2, R_2_r_2, delta_V_r, N_r ] = mirror_nonuniform(L_1_r_2, origin_1_r_2, mesh_1_r_2, cord_1, d,Lx,refinements); % Create discretization for mirror the refined parts of the cuboid
    [ R_2_r, delta_V_2_r, N_2_r ] = discretize_cuboid(L_2_r, origin_2_r, mesh_2_r, cord_2); % Create discretization for coarse part of cuboid #1
    [ R_2_c, delta_V_2_c, N_2_c ] = discretize_cuboid(L_2_c, origin_2_c, mesh_2_c, cord_2); % Create discretization for coarse part of cuboid #2
    r = [R_1_c;R_1_r; R_1_r_2; R_2_r_2; R_2_r; R_2_c]; % Combine discretization in one matrix
    N = N_1_c + N_1_r + N_r + N_2_r+ N_2_c; % Total number of subvolumes 
    N_1 = N_1_c + N_1_r + N_r/2;
    delta_V_vector = [delta_V_1_c, delta_V_1_r, delta_V_r, delta_V_2_r, delta_V_2_c]; % Vector of volume of each subvolume [m^3]
    
end   

ind_bulk = [1, N_1+1]; % Bulk object start index
L_sub = delta_V_vector.^(1/3); % Vector of length of a side of each cubic subvolumes

%%%%%%%%%%%%%
% File name %
%%%%%%%%%%%%%

back = cd;
saveDir = 'Discretizations/';

% File name for saved discretizations
file_name_saved = [description '_Lx' num2str((1e9)*Lx) 'nm_Ly' num2str((1e9)*Ly) 'nm_Lz' num2str((1e9)*Lz) 'nm_d' num2str((1e9)*d) 'nm_N'  num2str(N)]; % File name where results will be saved
   
%%%%%%%%%%%%%%%%%%%%%%%
% Plot discretization %
%%%%%%%%%%%%%%%%%%%%%%%

% Figure title
 title_string = {['Location of each subvolume for N = ', num2str(N) ' total subvolumes'],['L_x = ' num2str((1e6)*Lx) ' \mum, L_y = ' num2str((1e6)*Ly) ' \mum, L_z = ' num2str((1e9)*Lz) ' nm, d = ' num2str((1e9)*d) ' nm']};
  
% Visualize discretized lattice
FIG_discretization = figure(1);
plot3(r(:,1), r(:,2), r(:,3), 'x')
if show_titles == 1
    title(title_string, 'fontsize', 14)
end
xlabel('x-axis [m]')
ylabel('y-axis [m]')
zlabel('z-axis [m]')
if show_axes == 0
    grid off
    axis off
end
set(gca, 'Fontsize', 20)
set(gca,'DataAspectRatio',[1 1 1])
view(3)
grid on

% Plot voxel image of discretization
FIG_voxel = figure(2);
if refinements == 0  
    [vert, fac] = voxel_image( r, L_sub(1), [], [], [], [], 'on', [], [] );
end
if refinements == 1  
    [vert, fac] = voxel_image( R_1_c, L_sub(1), [], [], [], [], 'on', [], [] ); %L_sub(1) vox_sz - 1 x 3 vector with voxel size - if vox_sz is a scalar, all edges will have the same length
    [vert, fac] = voxel_image( R_1_r, L_sub(N_1_c+1), [], [], [], [], 'on', [], [] ); %L_sub(1) vox_sz - 1 x 3 vector with voxel size - if vox_sz is a scalar, all edges will have the same length
    [vert, fac] = voxel_image( R_2_r, L_sub(N_1+1), [], [], [], [], 'on', [], [] ); %L_sub(1) vox_sz - 1 x 3 vector with voxel size - if vox_sz is a scalar, all edges will have the same length
    [vert, fac] = voxel_image( R_2_c, L_sub(1), [], [], [], [], 'on', [], [] ); %L_sub(1) vox_sz - 1 x 3 vector with voxel size - if vox_sz is a scalar, all edges will have the same length
end

if refinements == 2  
    [vert, fac] = voxel_image( R_1_c, L_sub(1), [], [], [], [], 'on', [], [] ); %L_sub(1) vox_sz - 1 x 3 vector with voxel size - if vox_sz is a scalar, all edges will have the same length
    [vert, fac] = voxel_image( R_1_r, L_sub(N_1_c+1), [], [], [], [], 'on', [], [] ); %L_sub(1) vox_sz - 1 x 3 vector with voxel size - if vox_sz is a scalar, all edges will have the same length
    [vert, fac] = voxel_image( R_1_r_2, L_sub(N_1+1), [], [], [], [], 'on', [], [] ); %L_sub(1) vox_sz - 1 x 3 vector with voxel size - if vox_sz is a scalar, all edges will have the same length
    [vert, fac] = voxel_image( R_2_r_2, L_sub(N_1+1), [], [], [], [], 'on', [], [] ); %L_sub(1) vox_sz - 1 x 3 vector with voxel size - if vox_sz is a scalar, all edges will have the same length
    [vert, fac] = voxel_image( R_2_r, L_sub(N_1_c+1), [], [], [], [], 'on', [], [] ); %L_sub(1) vox_sz - 1 x 3 vector with voxel size - if vox_sz is a scalar, all edges will have the same length
    [vert, fac] = voxel_image( R_2_c, L_sub(1), [], [], [], [], 'on', [], [] ); %L_sub(1) vox_sz - 1 x 3 vector with voxel size - if vox_sz is a scalar, all edges will have the same length
end

if show_titles == 1
    title(title_string, 'fontsize', 14)
end
xlabel('x-axis (m)');
ylabel('y-axis (m)');
zlabel('z-axis (m)');
if show_axes == 0
    grid off
    axis off
end
set(gca, 'fontsize', 20)
view(35,20)
%view(5,20)
%view(2)


%%%%%%%%%%%%%%%%
% Save figures %
%%%%%%%%%%%%%%%%

if save_fig == 1
    % Center of subvolume discretization
    fig_path_1 = [saveDir 'figure_voxel/' file_name_saved '_discretization.fig'];
    saveas(FIG_discretization, fig_path_1)
    clear FIG_discretization % Remove previous plot handles

    
    % Voxel discretization
    fig_path_2 = [saveDir 'figure_voxel/' file_name_saved '_voxel.fig'];
    saveas(FIG_voxel, fig_path_2)
    clear FIG_voxel % Remove previous plot handles
    
end


%%%%%%%%%%%%%%%%%%%%%%%%
% Save discretizations %
%%%%%%%%%%%%%%%%%%%%%%%%

% Save discretization matrix to an Excel file
if save_Excel == 1
    disc_path_1 = [saveDir 'Excel_files/' file_name_saved '_discretization.xlsx'];
    writematrix(r, disc_path_1)
end

if save_txt == 1
% Save discretization matrix to a .txt file
disc_path_2 = [saveDir file_name_saved '_discretization.txt'];
writematrix(r, disc_path_2,'Delimiter','tab');

disc_path_3 = [saveDir file_name_saved '_delta_V_vector.txt'];
writematrix(delta_V_vector', disc_path_3,'Delimiter',',');
end


