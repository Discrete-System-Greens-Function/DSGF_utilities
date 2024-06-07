%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lindsay Walter and Livia Correa
% Thin film discretization
% Updated 03/12/23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear, clc, close all
format long
profile on
disp(['Running MATLAB script ' mfilename])

% This code creates discretizations for a system of two thin films
% separated by a gap.  Results are output in either .txt format or in an
% Excel file based on user input.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%******************************USER INPUTS********************************%

%%%%%%%%%%%%%%%%%%%%%%
% Set results export %
%%%%%%%%%%%%%%%%%%%%%%

back = cd;
saveDir = 'Discretizations/'; %fullfile(back, '../../../library/discretizations/thin-film/');

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


%%%%%%%%%%%%%%%%%%%%
% Dimension inputs %
%%%%%%%%%%%%%%%%%%%%

d = 100e-9;  % Distance between the two films [m]
Lx = 1000e-9;
Ly = 1000e-9;
Lz = 1000e-9;
refinements = 2;
slope = 'n';
rotate = 'n';
translate = 'n';


%%%%%%%%%%%%%%%%%%%%
% for t=20nm %
%%%%%%%%%%%%%%%%%%%%
%Lx_ref = 200e-9; % if refinement = 1 
Lx_ref = 180e-9; % if refinement = 2
Lx_ref_2 = 20e-9; % if refinement = 2 


%{
%%%%%%%%%%%%%%%%%%%%
% for t=50nm %
%%%%%%%%%%%%%%%%%%%%
%Lx_ref = 150e-9; % if refinement = 1 or 200e-9
Lx_ref = 150e-9; % if refinement =  2 
Lx_ref_2 = 50e-9; % if refinement = 2 
%}

%{
%%%%%%%%%%%%%%%%%%%%
% for t=120nm %
%%%%%%%%%%%%%%%%%%%%
%Lx_ref = 150e-9;
Lx_ref = 60e-9; % if refinement =  2 
Lx_ref_2 = 90e-9;% if refinement = 2  
%}


if refinements == 0
    % Thin film #1
    Lx_1 = Lx;         % Length of film #1 
    Ly_1 = Ly;         % Width of film #1
    Lz_1 = Lz;         % Thickness of film #1 
    origin_1 = [0,0,0];    % Point of back, left, bottom of cube [x,y,z]
    cord_1 = 'min';          % String of coordinate along which direction mesh is to be refined 
    mesh_1 = 1;            % Number of subvolumes across the dimension given in the input 'cord' 

    % Thin film #2
    Lx_2 = Lx_1;              % Length of film #2
    Ly_2 = Ly_1;              % Width of film #2
    Lz_2 = Lz_1;              % Thickness of film #2
    origin_2 = [Lx_1+d,0,0];  % Point of back, left, bottom of cube [x,y,z]
    cord_2 = 'min';             % String of coordinate along which direction mesh is to be refined 
    mesh_2 = mesh_1;               % Number of subvolumes across the dimension given in the input 'cord'   
end

if refinements == 1
    % Thin film #1
    Lx_1_c = Lx-Lx_ref;         % Length of film #1 
    Lx_1_r = Lx_ref;         % Length of film #1 
    Ly_1 = Ly;         % Width of film #1
    Lz_1 = Lz;         % Thickness of film #1 
    origin_1 = [0,0,0];    % Point of back, left, bottom of cube [x,y,z]
    origin_1_r = [Lx_1_c,0,0];    % Point of back, left, bottom of cube [x,y,z]
    cord_1 = 'min';          % String of coordinate along which direction mesh is to be refined 
    mesh_1_c = 1;            % Number of subvolumes across the dimension given in the input 'cord' 
    mesh_1_r = 2;            % Number of subvolumes across the dimension given in the input 'cord' 

    % Thin film #2
    Lx_2_c = Lx_1_c;              % Length of film #2
    Lx_2_r = Lx_1_r;              % Length of film #2
    Ly_2 = Ly_1;              % Width of film #2
    Lz_2 = Lz_1;              % Thickness of film #2
    origin_2_r = [d+Lx_1_c+Lx_1_r,0,0];  % Point of back, left, bottom of cube [x,y,z
    origin_2_c = [d+Lx_1_c+Lx_1_r+Lx_2_r,0,0];  % Point of back, left, bottom of cube [x,y,z]
    cord_2 = 'min';             % String of coordinate along which direction mesh is to be refined 
    mesh_2_c = mesh_1_c;               % Number of subvolumes across the dimension given in the input 'cord' 
    mesh_2_r = mesh_1_r;               % Number of subvolumes across the dimension given in the input 'cord'
end

if refinements == 2
    % Thin film #1
    Lx_1_c = Lx-Lx_ref-Lx_ref_2;         % Length of film #1 
    Lx_1_r = Lx_ref;         % Length of film #1 
    Lx_1_r_2 = Lx_ref_2;         % Length of film #1 
    Ly_1 = Ly;         % Width of film #1
    Lz_1 = Lz;         % Thickness of film #1 
    origin_1 = [0,0,0];    % Point of back, left, bottom of cube [x,y,z]
    origin_1_r = [Lx_1_c,0,0];    % Point of back, left, bottom of cube [x,y,z]
    origin_1_r_2 = [Lx_1_c+Lx_1_r,0,0];    % Point of back, left, bottom of cube [x,y,z]
    cord_1 = 'min';          % String of coordinate along which direction mesh is to be refined 
    mesh_1_c = 1;            % Number of subvolumes across the dimension given in the input 'cord' 
    mesh_1_r = 2;            % Number of subvolumes across the dimension given in the input 'cord' 
    mesh_1_r_2 = 4;            % Number of subvolumes across the dimension given in the input 'cord' 

    % Thin film #2
    Lx_2_c = Lx_1_c;              % Length of film #2
    Lx_2_r = Lx_1_r;              % Length of film #2
    Lx_2_r_2 = Lx_1_r_2;              % Length of film #2
    Ly_2 = Ly_1;              % Width of film #2
    Lz_2 = Lz_1;              % Thickness of film #2
    origin_2_r_2 = [d+Lx_1_c+Lx_1_r+Lx_1_r_2,0,0];    % Point of back, left, bottom of cube [x,y,z]
    origin_2_r = [d+Lx_1_c+Lx_1_r+Lx_1_r_2+Lx_2_r_2,0,0];  % Point of back, left, bottom of cube [x,y,z
    origin_2_c = [d+Lx_1_c+Lx_1_r+Lx_2_r+Lx_1_r_2+Lx_2_r_2,0,0];  % Point of back, left, bottom of cube [x,y,z]
    cord_2 = 'min';             % String of coordinate along which direction mesh is to be refined 
    mesh_2_c = mesh_1_c;               % Number of subvolumes across the dimension given in the input 'cord' 
    mesh_2_r = mesh_1_r;               % Number of subvolumes across the dimension given in the input 'cord'
    mesh_2_r_2 = mesh_1_r_2;               % Number of subvolumes across the dimension given in the input 'cord'
end

%****************************END USER INPUTS******************************%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%
% Create discretization %
%%%%%%%%%%%%%%%%%%%%%%%%%

if refinements == 0  
    % Combine inputs
    L_1 = [Lx_1, Ly_1, Lz_1]; % Vector containing length of each side of film #1 [m]
    L_2 = [Lx_1, Ly_1, Lz_1]; % Vector containing length of each side of film #1 [m]
    
    [ R_1, delta_V_1, N_1 ] = rectangular_prism_discretization(L_1, origin_1, mesh_1, cord_1);  % Create discretization for flat surface #1
    [ R_2, delta_V_2, N_2 ] = rectangular_prism_discretization(L_2, origin_2, mesh_2, cord_2); % Create discretization for flat surface #2
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
    
    [ R_1_c, delta_V_1_c, N_1_c ] = rectangular_prism_discretization(L_1_c, origin_1, mesh_1_c, cord_1); % Create discretization for coarse part of flat surface #1
    [ R_1_r, R_2_r, delta_V_r, N_r ] = edge_mirror_discretization(L_1_r, origin_1_r, mesh_1_r, cord_1, d, slope,Lx,refinements); % Create discretization for mirror the refined parts of the flat surfaces
    [ R_2_c, delta_V_2_c, N_2_c ] = rectangular_prism_discretization(L_2_c, origin_2_c, mesh_2_c, cord_2); % Create discretization for coarse part of flat surface #2
    r = [R_1_c;R_1_r; R_2_r; R_2_c]; % Combine discretization in one matrix
    N = N_1_c + N_r + N_2_c; % Total number of subvolumes 
    N_1 = N_1_c + N_r/2;
    delta_V_vector = [delta_V_1_c, delta_V_r, delta_V_2_c]; % Vector of volume of each subvolume [m^3]
    if slope == 'y'
    top = max(R_1_r(:,3));
    bot = min(R_2_r(:,3));

    Lsub_ref=min(Lx_1_r,Lz_1)/mesh_1_r;
    d_top = min(R_2_r(find(R_2_r(:,3) == top),1)) - max(R_1_r(find(R_1_r(:,3) == top),1))-Lsub_ref;
    d_bot = min(R_2_r(find(R_2_r(:,3) == bot),1)) - max(R_1_r(find(R_1_r(:,3) == bot),1))-Lsub_ref;
    %explanation: we verify that the dimensions of the membranes match what
    %we want (voxel plot). The calculation uses the discretization (star
    %plot) and we need to remove half of the subvolume in each
    %side=Lsub_ref.
    end
end    

if refinements == 2
    % Combine inputs
    L_1_c = [Lx_1_c, Ly_1, Lz_1]; % Vector containing length of each side of film #1 [m]
    L_1_r = [Lx_1_r, Ly_1, Lz_1]; % Vector containing length of each side of film #1 [m]
    L_1_r_2 = [Lx_1_r_2, Ly_1, Lz_1]; % Vector containing length of each side of film #1 [m]
    L_2_r_2 = [Lx_2_r_2, Ly_2, Lz_2]; % Vector containing length of each side of film #1 [m]
    L_2_r = [Lx_2_r, Ly_2, Lz_2]; % Vector containing length of each side of film #1 [m]
    L_2_c = [Lx_2_c, Ly_2, Lz_2]; % Vector containing length of each side of film #1 [m]
    
    [ R_1_c, delta_V_1_c, N_1_c ] = rectangular_prism_discretization(L_1_c, origin_1, mesh_1_c, cord_1); % Create discretization for coarse part of flat surface #1
    [ R_1_r, delta_V_1_r, N_1_r ] = rectangular_prism_discretization(L_1_r, origin_1_r, mesh_1_r, cord_1); % Create discretization for coarse part of flat surface #1
    [ R_1_r_2, R_2_r_2, delta_V_r, N_r ] = edge_mirror_discretization(L_1_r_2, origin_1_r_2, mesh_1_r_2, cord_1, d, slope,Lx,refinements); % Create discretization for mirror the refined parts of the flat surfaces
    [ R_2_r, delta_V_2_r, N_2_r ] = rectangular_prism_discretization(L_2_r, origin_2_r, mesh_2_r, cord_2); % Create discretization for coarse part of flat surface #1
    [ R_2_c, delta_V_2_c, N_2_c ] = rectangular_prism_discretization(L_2_c, origin_2_c, mesh_2_c, cord_2); % Create discretization for coarse part of flat surface #2
    r = [R_1_c;R_1_r; R_1_r_2; R_2_r_2; R_2_r; R_2_c]; % Combine discretization in one matrix
    N = N_1_c + N_1_r + N_r + N_2_r+ N_2_c; % Total number of subvolumes 
    N_1 = N_1_c + N_1_r + N_r/2;
    delta_V_vector = [delta_V_1_c, delta_V_1_r, delta_V_r, delta_V_2_r, delta_V_2_c]; % Vector of volume of each subvolume [m^3]
    if slope == 'y'
    
    Lsub_ref=min(Lz_1)/mesh_1_r_2;    
    
    bot = min(R_2_r_2(:,3)); 
    d_bot = min(R_2_r_2(find(R_2_r_2(:,3) == bot),1)) - max(R_1_r_2(find(R_1_r_2(:,3) == bot),1))-Lsub_ref;
   
    top = max(R_1_r_2(:,3));
    d_top = min(R_2_r_2(find(R_2_r_2(:,3) == top),1)) - max(R_1_r_2(find(R_1_r_2(:,3) == top),1))-Lsub_ref;
    %{
    top = max(R_1_r(:,3));
    Lsub_ref_coarse = min(Lx_1_r,Lz_1)/mesh_1_r;    
    d_top = min(R_2_r(find(R_2_r(:,3) == top),1)) - max(R_1_r(find(R_1_r(:,3) == top),1))-Lsub_ref_coarse;
    %}
    end
end   

ind_bulk = [1, N_1+1]; % Bulk object start index
L_sub = delta_V_vector.^(1/3); % Vector of length of a side of each cubic subvolumes

if rotate == 'y'
% Linear transformation for angular rotation of the membranes
% Rotate coordinates by 45 deg clockwise: https://www.mathworks.com/matlabcentral/answers/1653165-how-can-i-rotate-a-matrix-45-degrees?s_tid=ans_recom_leaf
angle = 5;
membrane_1= [R_1(:,1), R_1(:,3)]; %[1, 0];
rot_cw = [cosd(angle), sind(angle); -sind(angle), cosd(angle)];
x_z_rot_1 = membrane_1* rot_cw;  
R_1_rotate = [x_z_rot_1(:,1), R_1(:,2), x_z_rot_1(:,2)];

membrane_2= [R_2(:,1), R_2(:,3)]; %[1, 0];
rot_ccw = [cosd(angle), -sind(angle); sind(angle), cosd(angle)];
x_z_rot_2 = membrane_2* rot_ccw;  
%R_2_rotate = [x_z_rot_2(:,1), R_2(:,2), x_z_rot_2(:,2)+5.5e-7]; %include a translation of 5.5e-7 for 15o
R_2_rotate = [x_z_rot_2(:,1), R_2(:,2), x_z_rot_2(:,2)+1.8e-7]; %include a translation of 1.8e-7 for 5o


% Modified discretization with linear transformation
r = [R_1_rotate; R_2_rotate]; % combine discretization 

end

if translate == 'y'
%translation of one of the membranes
R_2_translate = [R_2(:,1), R_2(:,2), R_2(:,2)+1.8e-7]; 
% Modified discretization with linear transformation
r = [R_1; R_2_translate]; % combine discretization 
end




%%%%%%%%%%%%%
% File name %
%%%%%%%%%%%%%

% File name for saved discretizations
description = '2_thin_films';  % Very short description of results
%file_name_saved = [description '_Lx' num2str((1e6)*Lx_1) 'um_Ly' num2str((1e6)*Ly_1) 'um_Lz' num2str((1e9)*Lz_1) 'nm_d' num2str((1e9)*d) 'nm_N'  num2str(N)]; % File name where results will be saved
if slope == 'n'
    file_name_saved = [description '_Lx' num2str((1e9)*Lx) 'nm_Ly' num2str((1e9)*Ly) 'nm_Lz' num2str((1e9)*Lz) 'nm_d' num2str((1e9)*d) 'nm_N'  num2str(N)]; % File name where results will be saved
elseif slope == 'y'
    file_name_saved = [description '_Lx' num2str((1e9)*Lx) 'nm_Ly' num2str((1e9)*Ly) 'nm_Lz' num2str((1e9)*Lz) 'nm_d_min' num2str((1e9)*d_bot) 'nm_d_max' num2str((1e9)*d_top) 'nm_N'  num2str(N)]; % File name where results will be saved
end    
%%%%%%%%%%%%%%%%%%%%%%%
% Plot discretization %
%%%%%%%%%%%%%%%%%%%%%%%

% Figure title
if slope == 'n' 
 title_string = {['Location of each subvolume for N = ', num2str(N) ' total subvolumes'],['L_x = ' num2str((1e6)*Lx) ' \mum, L_y = ' num2str((1e6)*Ly) ' \mum, L_z = ' num2str((1e9)*Lz) ' nm, d = ' num2str((1e9)*d) ' nm']};
elseif slope == 'y'
 title_string = {['Location of each subvolume for N = ', num2str(N) ' total subvolumes'],['L_x = ' num2str((1e6)*Lx) ' \mum, L_y = ' num2str((1e6)*Ly) ' \mum, L_z = ' num2str((1e9)*Lz) ' nm, d_{min} = ' num2str((1e9)*d_bot) ' nm, d_{max} = ' num2str((1e9)*d_top) ' nm']};
end    
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


