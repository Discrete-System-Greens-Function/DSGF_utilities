%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Livia Correa and Lindsay Walter
% Discretization for two user-defined thermal sources 
% Updated 06/27/24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear Workspace and close all figures
clear, clc, close all

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Import discretizations %
%%%%%%%%%%%%%%%%%%%%%%%%%%

discDir = "single_sources/"; % Directory where discretization 1 is stored

discretization = ["cube_125.xlsx","cube_125.xlsx"]; % File name of discretization 
L_char = [500.e-9; 500.e-9];
d = 500.e-9; 

% Define origin depending of the source shape. Edge-to-edge gap spacing is used.

% Default:
%origin = [0, 0, 0; d + L_char(1) + L_char(2) , 0, 0]; 

% Uncomment if second source is a cube:
origin = [0, 0, 0; d + L_char(1)/2 + L_char(2)/2 , 0, 0]; 

%%%%%%%%%%%%%%%%%%%%%%
% Set results export %
%%%%%%%%%%%%%%%%%%%%%%

back = cd;
saveDir = 'Discretizations/'; 

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

% reads all the discretizations for all the bulk objects
%
%	Inputs:
%		discretization - vector of the all the bulk objects and their discretizations
%		L_char - vector of the characteristic lengths of the objects
%		origin - vector of the origins of all the bulk objects
%
%	Outputs:
%		N_each_object
%		volume
%		r_each_object
%		ind_bulk
%		delta_V_each_object
%		L_sub_each_object

N_bulk = length(discretization); % Number of bulk objects

% Determine file structure of bulk object discretization and extract discretization information.
N_each_object = zeros(N_bulk,1);       % Preallocate
volume = zeros(N_bulk,1);              % Preallocate
r_each_object = cell(N_bulk,1);        % Preallocate
ind_bulk = ones(N_bulk,1);             % Preallocate
delta_V_each_object = cell(N_bulk,1);  % Preallocate
L_sub_each_object = cell(N_bulk,1);    % Preallocate
geometry = cell(N_bulk,1);        % Preallocate
geometry_N = cell(N_bulk,1);        % Preallocate
delta_V = zeros(N_bulk,1);              % Preallocate
L_sub = zeros(N_bulk,1);              % Preallocate
for ii = 1:2
       
    %discFile = string(append(discretization{ii},file_format));
    discFile = discretization{ii};
    geometry_N{ii} = extractBefore(discFile, '.');     % Geometry of bulk object
    geometry{ii} = extractBefore(discFile, '_');     % Geometry of bulk object
    
    % Import unscaled discretization of each object
    r_each_object{ii} = readmatrix(append(discDir, discFile));

    % Number of subvolumes in each bulk object
    [N_each_object(ii),~] = size(r_each_object{ii});

    
    % Subvolume size for each object (uniform discretization)
    % Scale sample discretizations based on the geometryswitch (geometry)
    
	if strcmp("sphere",geometry{ii})
		volume(ii) = (4/3)*pi*(L_char(ii)^3);  % Volume of sphere [m^3]
        delta_V(ii) = (volume(ii)/N_each_object(ii)); % Volume of subvolumes in each object (uniform discretization)
        L_sub(ii) = delta_V(ii)^(1/3); % Length of side of a cubic subvolume in each object (uniform discretization)
        delta_V_each_object{ii} = ones(N_each_object(ii), 1).*delta_V(ii); % Volume of subvolumes in each object (uniform discretization)
        L_sub_each_object{ii} = ones(N_each_object(ii), 1).*L_sub(ii);  % Length of side of a cubic subvolume in each object (uniform discretization)
	elseif strcmp("cube",geometry{ii})
		volume(ii) = (L_char(ii)^3);           % Volume of cube [m^3]
        delta_V(ii) = (volume(ii)/N_each_object(ii)); % Volume of subvolumes in each object (uniform discretization)
        L_sub(ii) = delta_V(ii)^(1/3);
        delta_V_each_object{ii} = ones(N_each_object(ii), 1).*delta_V(ii); % Volume of subvolumes in each object (uniform discretization)
        L_sub_each_object{ii} = ones(N_each_object(ii), 1).*L_sub(ii);  % Length of side of a cubic subvolume in each object (uniform discretization)
	elseif strcmp("dipole",geometry{ii})
		volume(ii) = (4/3)*pi*(L_char.^3);  % Volume of spherical dipole [m^3]
        delta_V(ii) = (volume(ii)/N_each_object(ii)); % Volume of subvolumes in each object (uniform discretization)
        L_sub(ii) = delta_V(ii)^(1/3); % Length of side of a cubic subvolume in each object (uniform discretization)
        delta_V_each_object{ii} = ones(N_each_object, 1).*delta_V(ii); % Volume of subvolumes in each object (uniform discretization)
        L_sub_each_object{ii} = ones(N_each_object(ii), 1).*L_sub(ii);  % Length of side of a cubic subvolume in each object (uniform discretization) 
    else % general case (uniform discretization)
        [L_sub(ii), delta_V(ii)] = calculate_Lsub_uniform(r_each_object{ii}); % Calculates the subvolume side length (i.e., cubic lattice discretization length scale) for a uniform cubic lattice discretization
        delta_V_each_object{ii} = ones(N_each_object(ii), 1).*delta_V(ii); % Volume of subvolumes in each object (uniform discretization)
        L_sub_each_object{ii} = ones(N_each_object(ii), 1).*L_sub(ii);     % Length of side of a cubic subvolume in each object (uniform discretization)
	end % End switch-case through geometries
    
    if ii ~= 1
	    	ind_bulk(ii) = sum(N_each_object(1:ii-1))+1;
    end
    
    % Scale discretization
    r_each_object{ii} = L_sub_each_object{ii}.*r_each_object{ii};

	% Move the center-of-mass of each object to the origin [0,0,0]
	r_each_object{ii} = center_of_mass(r_each_object{ii});

    % Move each discretization to its user-specified origin
	r_each_object{ii} = r_each_object{ii} + repmat(origin(ii,:), N_each_object(ii), 1);
        
end

% Discretized lattice including subvolumes of all objects in one matrix (N x 3 matrix)
r = cell2mat(r_each_object);

% Total number of subvolumes
[N,~] = size(r);

% Subvolume size for all N subvolumes (N x 1 vectors)
delta_V_vector = cell2mat(delta_V_each_object); % Volume of subvolumes for all N subvolumes
L_sub_vector = cell2mat(L_sub_each_object);     % Length of side of a cubic subvolume for all N subvolumes

%%%%%%%%%%%%%
% File name %
%%%%%%%%%%%%%

file_name_saved = [geometry_N{1} '_Lchar_' num2str(L_char(1)) '_d_' num2str(d) '_' geometry_N{2} '_Lchar_' num2str(L_char(2))]; % File name where results will be saved

%%%%%%%%%%%%%%%%%%%%%%%
% Plot discretization %
%%%%%%%%%%%%%%%%%%%%%%%

title_string = {'N =' num2str(N)}; % Figure title
 
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot cubic lattice figure %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FIG_voxel = figure(2);
for ii = 1:2
[vert, fac] = voxel_image( r_each_object{ii}, L_sub(ii), [], [], [], [], 'on', [], [] );
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

% Save discretization and delta_V_vector as .txt files
if save_txt == 1
    disc_path_2 = [saveDir file_name_saved '_discretization.txt'];
    writematrix(r, disc_path_2,'Delimiter','tab');

    disc_path_3 = [saveDir file_name_saved '_delta_V_vector.txt'];
    writematrix(delta_V_vector, disc_path_3,'Delimiter',',');
end











