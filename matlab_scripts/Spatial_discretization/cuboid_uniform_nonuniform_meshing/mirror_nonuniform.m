function [ R_1, R_2, delta_V, N_r ] = mirror_nonuniform( L, origin, mesh, cord, d, Lx, refinements )

N = zeros(1,3); % Preallocate vector of number of discretizations for each side
L_min = L(3);%min(L);              % Find length of smallest side
delta_L = L_min/mesh;        % Size of discretized length
N = round(L./delta_L);       % Vector of number of discretizations for each side

% Cut cube into subvolumes
x1 = linspace(origin(1), origin(1)+L(1), N(1)+1);  % Location of discretized edges along x direction
y1 = linspace(origin(2), origin(2)+L(2), N(2)+1);  % Location of discretized edges along y direction
z1 = linspace(origin(3), origin(3)+L(3), N(3)+1);  % Location of discretized edges along z direction

    % Make discretized lattice for rectangular prism
    %R_1 = zeros(N_disc, 3);  % Preallocate location vector matrix
    n = 1;                 % Set number of iterations to 1
    for k = 1:N(3)         % Loop through discretized z axis
        for j = 1:N(2)     % Loop through discretized y axis
            for h = 1:N(1) % Loop through discretized x axis
                x = x1(h) + delta_L/2; % Location of discretized volume centers along x direction
                y = y1(j) + delta_L/2; % Location of discretized volume centers along y direction
                z = z1(k) + delta_L/2; % Location of discretized volume centers along z direction
                R_1(n,1) = x; % Populate location matrix
                R_1(n,2) = y; % Populate location matrix
                R_1(n,3) = z; % Populate location matrix
                n = n+1;    % Increase count of loop iteration
            end
        end
    end 

R_2 = [-R_1(:,1)+d+Lx*2,R_1(:,2),R_1(:,3)];
%%

N_r = length(R_1)+ length(R_2); % Total number of subvolumes in rectangular prism
delta_V = (delta_L^3)*ones(1,N_r); % Discretized volume (volume of a single subvolume)




