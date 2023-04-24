function [ R_1, R_2, delta_V, N_r ] = edge_mirror_discretization( L, origin, mesh, cord, d, slope, Lx, refinements )

N = zeros(1,3); % Preallocate vector of number of discretizations for each side
L_min = L(3);%min(L);              % Find length of smallest side
delta_L = L_min/mesh;        % Size of discretized length
N = round(L./delta_L);       % Vector of number of discretizations for each side

% Cut cube into subvolumes
x1 = linspace(origin(1), origin(1)+L(1), N(1)+1);  % Location of discretized edges along x direction
y1 = linspace(origin(2), origin(2)+L(2), N(2)+1);  % Location of discretized edges along y direction
z1 = linspace(origin(3), origin(3)+L(3), N(3)+1);  % Location of discretized edges along z direction

if slope == 'n'

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
end

if slope == 'y'
   if refinements == 1
    n=1;
    subvol_ref = L(1)/delta_L
    for k = 1:N(3)
    for j = 1:N(2)     % Loop through discretized y axis
        for i=1:subvol_ref -(k-1)
            x = x1(1) + (2*i-1)*delta_L/2; % Location of discretized volume centers along x direction delta_L/2
            y = y1(j)+ delta_L/2; % Location of discretized volume centers along y direction
            z = z1(k)+ delta_L/2; % Location of discretized volume centers along z direction
            R_1(n,1) = x; % Populate location matrix
            R_1(n,2) = y; % Populate location matrix
            R_1(n,3) = z; % Populate location matrix
            %-(k+1)
            n=n+1;
        end
    end
    end
    
    
        % this condition is used to fill an extra row with subvolumes, comment if not used 
        n_extra = n-1;
        %{
        if k== N(3)-1
        for j = 1:N(2)     % Loop through discretized y axis
        for i=1:2%subvol_ref -(k-1)
            x = R_1(n_extra,1) + i*delta_L; % Location of discretized volume centers along x direction delta_L/2
            y = y1(j)+ delta_L/2; % Location of discretized volume centers along y direction
            z = z1(k)+ delta_L/2; % Location of discretized volume centers along z direction
            R_1(n,1) = x; % Populate location matrix
            R_1(n,2) = y; % Populate location matrix
            R_1(n,3) = z; % Populate location matrix
            %-(k+1)
            n=n+1;
        end
        end
        end
        %}
    
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Lz = 50nm, 4 layers of subvolumes
        % Lz = 120nm, 8 layers of subvolumes
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if k== N(3)
        for j = 1:N(2)     % Loop through discretized y axis
        for i=1:1%subvol_ref -(k-1)
            x = R_1(n_extra,1) + i*delta_L; % Location of discretized volume centers along x direction delta_L/2
            y = y1(j)+ delta_L/2; % Location of discretized volume centers along y direction
            z = z1(k)+ delta_L/2; % Location of discretized volume centers along z direction
            R_1(n,1) = x; % Populate location matrix
            R_1(n,2) = y; % Populate location matrix
            R_1(n,3) = z; % Populate location matrix
            %-(k+1)
            n=n+1;
        end
        end
        end
        
        
       
       
    end % refinements
    
    if refinements == 2
    n=1;
    subvol_ref = L(1)/delta_L
        
    %general case    
    for k = 1:N(3)
    for j = 1:N(2)     % Loop through discretized y axis
        for i=1:subvol_ref -(k-1)
            x = x1(1) + (2*i-1)*delta_L/2; % Location of discretized volume centers along x direction delta_L/2
            y = y1(j)+ delta_L/2; % Location of discretized volume centers along y direction
            z = z1(k)+ delta_L/2; % Location of discretized volume centers along z direction
            R_1(n,1) = x; % Populate location matrix
            R_1(n,2) = y; % Populate location matrix
            R_1(n,3) = z; % Populate location matrix
            %-(k+1)
            n=n+1;
        end    
    end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Lz = 120nm, 8 layers of subvolumes
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if k == 9 | k == 10 
        for j = 1:N(2)     % Loop through discretized y axis
        for i=1:2%subvol_ref -(k-1)
            x =  x1(1) + (2*i-1)*delta_L/2; % Location of discretized volume centers along x direction delta_L/2
            y = y1(j)+ delta_L/2; % Location of discretized volume centers along y direction
            z = z1(k)+ delta_L/2; % Location of discretized volume centers along z direction
            R_1(n,1) = x; % Populate location matrix
            R_1(n,2) = y; % Populate location matrix
            R_1(n,3) = z; % Populate location matrix
            %-(k+1)
            n=n+1;
        end
        end
        end
        
        if k == 11 | k == 12
        for j = 1:N(2)     % Loop through discretized y axis
        for i=1:1%subvol_ref -(k-1)
            x =  x1(1) + (2*i-1)*delta_L/2; % Location of discretized volume centers along x direction delta_L/2
            y = y1(j)+ delta_L/2; % Location of discretized volume centers along y direction
            z = z1(k)+ delta_L/2; % Location of discretized volume centers along z direction
            R_1(n,1) = x; % Populate location matrix
            R_1(n,2) = y; % Populate location matrix
            R_1(n,3) = z; % Populate location matrix
            %-(k+1)
            n=n+1;
        end
        end
        end
        
        
    
    %{
    % This modification is used for t =120nm, mesh_1_c = 1, mesh_1_r = 4, mesh_1_r_2 = 8;  
    
    if k== N(3)
        for j = 1:N(2)     % Loop through discretized y axis
        for i=1:subvol_ref -(k-1)
            x = R_1(n_extra,1) + delta_L; % Location of discretized volume centers along x direction delta_L/2
            y = y1(j)+ delta_L/2; % Location of discretized volume centers along y direction
            z = z1(k)+ delta_L/2; % Location of discretized volume centers along z direction
            R_1(n,1) = x; % Populate location matrix
            R_1(n,2) = y; % Populate location matrix
            R_1(n,3) = z; % Populate location matrix
            %-(k+1)
            n=n+1;
        end
        end
   end
    %}
    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Lz = 50nm, N=15680, 8 layers of subvolumes
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for k=1:2
    for j = 1:N(2)     % Loop through discretized y axis
        for i=1:subvol_ref 
            x = x1(1) + (2*i-1)*delta_L/2; % Location of discretized volume centers along x direction delta_L/2
            y = y1(j)+ delta_L/2; % Location of discretized volume centers along y direction
            z = z1(k)+ delta_L/2; % Location of discretized volume centers along z direction
            R_1(n,1) = x; % Populate location matrix
            R_1(n,2) = y; % Populate location matrix
            R_1(n,3) = z; % Populate location matrix
            %-(k+1)
            n=n+1;
        end
    end
    end
    
    for k=3:4
    for j = 1:N(2)     % Loop through discretized y axis
        for i=1:subvol_ref -2
            x = x1(1) + (2*i-1)*delta_L/2; % Location of discretized volume centers along x direction delta_L/2
            y = y1(j)+ delta_L/2; % Location of discretized volume centers along y direction
            z = z1(k)+ delta_L/2; % Location of discretized volume centers along z direction
            R_1(n,1) = x; % Populate location matrix
            R_1(n,2) = y; % Populate location matrix
            R_1(n,3) = z; % Populate location matrix
            %-(k+1)
            n=n+1;
        end
    end
    end
    
     for k = 5:8
    for j = 1:N(2)     % Loop through discretized y axis
        for i=1:subvol_ref -4
            x = x1(1) + (2*i-1)*delta_L/2; % Location of discretized volume centers along x direction delta_L/2
            y = y1(j)+ delta_L/2; % Location of discretized volume centers along y direction
            z = z1(k)+ delta_L/2; % Location of discretized volume centers along z direction
            R_1(n,1) = x; % Populate location matrix
            R_1(n,2) = y; % Populate location matrix
            R_1(n,3) = z; % Populate location matrix
            %-(k+1)
            n=n+1;
        end
    end
     end
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % end  Lz = 50nm, N=15680, 8 layers of subvolumes
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %}
     
    %{
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Lz = 120nm, N=, 8 layers of subvolumes
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for k=1:2
    for j = 1:N(2)     % Loop through discretized y axis
        for i=1:subvol_ref 
            x = x1(1) + (2*i-1)*delta_L/2; % Location of discretized volume centers along x direction delta_L/2
            y = y1(j)+ delta_L/2; % Location of discretized volume centers along y direction
            z = z1(k)+ delta_L/2; % Location of discretized volume centers along z direction
            R_1(n,1) = x; % Populate location matrix
            R_1(n,2) = y; % Populate location matrix
            R_1(n,3) = z; % Populate location matrix
            %-(k+1)
            n=n+1;
        end
    end
    end
    
    for k=3:4
    for j = 1:N(2)     % Loop through discretized y axis
        for i=1:subvol_ref -2
            x = x1(1) + (2*i-1)*delta_L/2; % Location of discretized volume centers along x direction delta_L/2
            y = y1(j)+ delta_L/2; % Location of discretized volume centers along y direction
            z = z1(k)+ delta_L/2; % Location of discretized volume centers along z direction
            R_1(n,1) = x; % Populate location matrix
            R_1(n,2) = y; % Populate location matrix
            R_1(n,3) = z; % Populate location matrix
            %-(k+1)
            n=n+1;
        end
    end
    end
    
    for k = 5:6%3:4
    for j = 1:N(2)     % Loop through discretized y axis
        for i=1:subvol_ref -2
            x = x1(1) + (2*i-1)*delta_L/2; % Location of discretized volume centers along x direction delta_L/2
            y = y1(j)+ delta_L/2; % Location of discretized volume centers along y direction
            z = z1(k)+ delta_L/2; % Location of discretized volume centers along z direction
            R_1(n,1) = x; % Populate location matrix
            R_1(n,2) = y; % Populate location matrix
            R_1(n,3) = z; % Populate location matrix
            %-(k+1)
            n=n+1;
        end
    end
    end
    
    for k = 7:8%3:4
    for j = 1:N(2)     % Loop through discretized y axis
        for i=1:subvol_ref -3
            x = x1(1) + (2*i-1)*delta_L/2; % Location of discretized volume centers along x direction delta_L/2
            y = y1(j)+ delta_L/2; % Location of discretized volume centers along y direction
            z = z1(k)+ delta_L/2; % Location of discretized volume centers along z direction
            R_1(n,1) = x; % Populate location matrix
            R_1(n,2) = y; % Populate location matrix
            R_1(n,3) = z; % Populate location matrix
            %-(k+1)
            n=n+1;
        end
    end
    end
    
    for k = 9:10%5:6
    for j = 1:N(2)     % Loop through discretized y axis
        for i=1:subvol_ref -4
            x = x1(1) + (2*i-1)*delta_L/2; % Location of discretized volume centers along x direction delta_L/2
            y = y1(j)+ delta_L/2; % Location of discretized volume centers along y direction
            z = z1(k)+ delta_L/2; % Location of discretized volume centers along z direction
            R_1(n,1) = x; % Populate location matrix
            R_1(n,2) = y; % Populate location matrix
            R_1(n,3) = z; % Populate location matrix
            %-(k+1)
            n=n+1;
        end
    end
    end
    
    for k = 11:12%5:6
    for j = 1:N(2)     % Loop through discretized y axis
        for i=1:1 %subvol_ref -4
            x = x1(1) + (2*i-1)*delta_L/2; % Location of discretized volume centers along x direction delta_L/2
            y = y1(j)+ delta_L/2; % Location of discretized volume centers along y direction
            z = z1(k)+ delta_L/2; % Location of discretized volume centers along z direction
            R_1(n,1) = x; % Populate location matrix
            R_1(n,2) = y; % Populate location matrix
            R_1(n,3) = z; % Populate location matrix
            %-(k+1)
            n=n+1;
        end
    end
    end
    
    
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % end  Lz = 120nm, N=, 8 layers of subvolumes
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %}
     
    
    
    %{
    for k = 7:8
    for j = 1:N(2)     % Loop through discretized y axis
        for i=1:2 %subvol_ref -4
            x = x1(1) + (2*i-1)*delta_L/2; % Location of discretized volume centers along x direction delta_L/2
            y = y1(j)+ delta_L/2; % Location of discretized volume centers along y direction
            z = z1(k)+ delta_L/2; % Location of discretized volume centers along z direction
            R_1(n,1) = x; % Populate location matrix
            R_1(n,2) = y; % Populate location matrix
            R_1(n,3) = z; % Populate location matrix
            %-(k+1)
            n=n+1;
        end
    end
    end
    %}
    
   
    
    
        
        %{
        % This modification is used for t =50nm, mesh_1_c = 1, mesh_1_r = 2, mesh_1_r_2 = 4; 
        n_extra = n-1;
        if k==N(3)
        for j = 1:N(2)     % Loop through discretized y axis
        for i=1:subvol_ref -(k-1)
            x = R_1(n_extra,1) + delta_L; % Location of discretized volume centers along x direction delta_L/2
            y = y1(j)+ delta_L/2; % Location of discretized volume centers along y direction
            z = z1(k)+ delta_L/2; % Location of discretized volume centers along z direction
            R_1(n,1) = x; % Populate location matrix
            R_1(n,2) = y; % Populate location matrix
            R_1(n,3) = z; % Populate location matrix
            %-(k+1)
            n=n+1;
        end
        end
        end
        %}
    
        %{
        for j = 1:N(2)     % Loop through discretized y axis
        for i=1:subvol_ref 
            x = R_1(n,1) + delta_L; % Location of discretized volume centers along x direction delta_L/2
            y = y1(j)+ delta_L/2; % Location of discretized volume centers along y direction
            z = z1(2)+ delta_L/2; % Location of discretized volume centers along z direction
            R_1(n,1) = x; % Populate location matrix
            R_1(n,2) = y; % Populate location matrix
            R_1(n,3) = z; % Populate location matrix
            %-(k+1)
            n=n+1;
        end
        end
        %}
    
        %{
        n_extra = n-1;
        if k== N(3)
        for j = 1:N(2)     % Loop through discretized y axis
        for i=1:1%subvol_ref -(k-1)
            x = R_1(n_extra,1) + i*delta_L; % Location of discretized volume centers along x direction delta_L/2
            y = y1(j)+ delta_L/2; % Location of discretized volume centers along y direction
            z = z1(k)+ delta_L/2; % Location of discretized volume centers along z direction
            R_1(n,1) = x; % Populate location matrix
            R_1(n,2) = y; % Populate location matrix
            R_1(n,3) = z; % Populate location matrix
            %-(k+1)
            n=n+1;
        end
        end
        end
        %}
    
    end %for k
    
    
    end % refinements
    
end %slope
    

R_2 = [-R_1(:,1)+d+Lx*2,R_1(:,2),R_1(:,3)];
%%

N_r = length(R_1)+ length(R_2); % Total number of subvolumes in rectangular prism
delta_V = (delta_L^3)*ones(1,N_r); % Discretized volume (volume of a single subvolume)




