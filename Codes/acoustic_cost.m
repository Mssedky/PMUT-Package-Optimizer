% Cost function evaluation for PMUT packaging for acoustic pressure 
% enhacement using a k-Wave acoustic simulation model
% Code by Mostafa Sedky
% August 1, 2026

function cost = acoustic_cost(lambda)
    %% ========== PARAMETERS ==========
    % PMUT parameters
    pmut_diameter = 0.7e-3;           % pMUT diameter (m)
    pmut_radius = pmut_diameter/2;
    f0 = 180e3;                     % Resonance frequency (Hz)
    v0 = 1e-5;                      % Peak particle velocity (m/s)
    base_amplitude = 106.77*1e-6;    % Base amplitude at 4 mm away from source with no horn (but stock short cylinder) (Pa)
    
    % Acoustic housing (horn) parameters
    cylinder_diameter = 8e-3;       % Outer cylinder diameter (m)
    cylinder_radius = cylinder_diameter/2;
    cylinder_height = 0.475e-3;      % Base packaging cylinder height (m) (should be 0.475 but can add 0.1 mm for printed housing)
    bottom_cylinder_radius = pmut_radius * 1.2;
    makeHorn = 1;
        
    % Generalized ellipsoidal hole shape parameters
    horn_shape.r_base = pmut_radius * 1.2 * 2;  % Base radius at PMUT

    % OPTIMIZATION PARAMETERS:
    horn_height = lambda(1);             % Height of acoustic housing (m)
    horn_shape.r_top = lambda(2);                % Top radius
    horn_shape.p = lambda(3);                        % Ellipsoid exponent p
    horn_shape.q = lambda(4);                        % Ellipsoid exponent q
    horn_shape.r = lambda(5);                        % Ellipsoid exponent r
    horn_shape.A_theta = lambda(6);       % Azimuthal modulation amplitude
    horn_shape.m_theta = lambda(7);         % Azimuthal mode number
    horn_shape.A_x = lambda(8);           % X-direction modulation
    horn_shape.A_y = lambda(9);           % Y-direction modulation
    horn_shape.A_z = lambda(10);           % Z-direction modulation
    horn_shape.m_z = lambda(11);             % Z-direction mode number
    
    % Medium properties
    c_air = 343;                    % Speed of sound in air (m/s)
    rho_air = 1.21;                 % Density of air (kg/m^3)
    c_horn = 2730;                  % Speed of sound in horn material (m/s) - e.g., aluminum
    rho_horn = 1190;                % Density of horn material (kg/m^3)
    
    % Simulation parameters
    z_up_target = 4e-3;             % Target height to measure pressue (m)
    ppw = 30;                       % Points per wavelength
    pml_size = 5;                   % PML Layer size 
    
    Nsteps = 5000;                  % Total number of time steps (2185)
    deltaT = 4e-9;                  % time step size (s)
    
    
    %% ========== COMPUTATIONAL GRID SETUP ==========
    lambda = c_air / f0;
    dx = lambda / ppw;
    dy = dx;
    dz = dx;

    % ---- User-defined physical grid size ----
    grid_length_x = 5.5e-3;     % Total width of domain in x [m]
    grid_length_y = 5.5e-3;     % Total width of domain in y [m]
    grid_length_z = 5e-3;     % Total height (z-direction) [m]
    
    % ---- Compute grid points ----
    Nx = round(grid_length_x / dx);
    Ny = round(grid_length_y / dy);
    Nz = round(grid_length_z / dz);
    
    % Make dimensions even for efficiency
    Nx = Nx + mod(Nx, 2);
    Ny = Ny + mod(Ny, 2);
    Nz = Nz + mod(Nz, 2);
   
    kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);
    
    % Time array
    kgrid.Nt = Nsteps;
    kgrid.dt = deltaT;
   
    %% ========== GENERATE ACOUSTIC HOUSING GEOMETRY ==========

    % Create coordinate grids
    x = kgrid.x_vec;
    y = kgrid.y_vec;
    z = kgrid.z_vec;
    [X_grid, Y_grid, Z_grid] = meshgrid(x, y, z);
    
    % Place pMUT right after PML at the bottom interior cell:
    source_z_index = pml_size + 1;    % first non-PML z-index (bottom interior)
    z_base = z(source_z_index) + cylinder_height;       % physical z coordinate of the base of the horn / pMUT plane
    
    % Initialize medium property matrices
    sound_speed = c_air * ones(Ny, Nx, Nz);
    density = rho_air * ones(Ny, Nx, Nz);
    
    % Generate horn profile
    n_z_horn = round(horn_height / dz);
    z_horn = linspace(0, horn_height, n_z_horn);
    n_theta = 360;
    theta = linspace(0, 2*pi, n_theta);
    
    % Compute base ellipsoidal radius profile
    r_curve = zeros(size(z_horn));
    for i = 1:length(z_horn)
        % Generalized ellipsoid interpolation
        t = z_horn(i) / horn_height;  % Normalized height [0,1]
        r_curve(i) = horn_shape.r_base * (1-t)^(1/horn_shape.p) + ...
                     horn_shape.r_top * t^(1/horn_shape.q);
    end
    
    % Apply modulations to create complex shape
    scale = @(idx) (z_horn(idx) / horn_height).^horn_shape.r;  % Normalized height
    kx = 2*pi / horn_height;
    ky = 2*pi / horn_height;
    kz = 2*pi / horn_height;
    
    horn_X = zeros(length(z_horn), n_theta);
    horn_Y = zeros(length(z_horn), n_theta);
    horn_Z = zeros(length(z_horn), n_theta);
    
    for i = 1:length(z_horn)
        mod_theta = 1 + (horn_shape.A_theta * scale(i)) * sin(horn_shape.m_theta * theta);
        mod_x = 1 + horn_shape.A_x * scale(i) * sin(kx * z_horn(i));
        mod_y = 1 + horn_shape.A_y * scale(i) * sin(ky * z_horn(i));
        
        z_mod = z_horn(i);
        r_i = r_curve(i) * mod_theta .* (1 + horn_shape.A_z * scale(i) * sin(horn_shape.m_z * kz * z_horn(i)));
        
        horn_X(i,:) = r_i .* cos(theta) * mod_x;
        horn_Y(i,:) = r_i .* sin(theta) * mod_y;
        horn_Z(i,:) = z_mod;
    end
    
    % Create binary mask for horn solid regions (these will be rigid boundaries)
    horn_mask = zeros(Ny, Nx, Nz);
    if makeHorn == 1
        for iz = 1:Nz
            z_current = z(iz);
            
            if z_current >= z_base && z_current <= (z_base + horn_height)
                % relative height inside horn (0 .. horn_height)
                z_rel = z_current - z_base;
        
                % Find corresponding horn profile
                [~, z_idx] = min(abs(z_horn - z_rel));
                r_hole = max(sqrt(horn_X(z_idx,:).^2 + horn_Y(z_idx,:).^2));
                
                for ix = 1:Nx
                    for iy = 1:Ny
                        x_curr = x(ix);
                        y_curr = y(iy);
                        r = sqrt(x_curr^2 + y_curr^2);
                        
                        % Inside outer cylinder but outside hole = solid horn
                        if r <= cylinder_radius && r > r_hole
                            horn_mask(iy, ix, iz) = 1;
                        end
                    end
                end
            end
        end
    end
    
    % --- Add the solid cylinder below horn ---
    cylinder_region = (Z_grid < z_base) & (Z_grid >= z_base - cylinder_height) & ...
                      (sqrt(X_grid.^2 + Y_grid.^2) >= bottom_cylinder_radius);
    horn_mask = horn_mask | cylinder_region;
    

    %% ========== SOURCE: pMUT VELOCITY PROFILE ==========
    % For a circular clamped membrane vibrating in (0,1) mode:
    % v(r,t) = v0 * J0(k*r) * sin(omega*t)
    
    % Membrane parameters
    k_membrane = 2.4048 / pmut_radius;  % Wavenumber for (0,1) mode
    
    % --- Find horn base centroid (use horn profile at z = 0) ---
    base_x_coords = horn_X(1,:);   % coordinates of hole boundary in x (m)
    base_y_coords = horn_Y(1,:);   % coordinates of hole boundary in y (m)
    
    % Compute centroid of the hole boundary (in physical coords)
    centroid_hole_x = mean(base_x_coords);
    centroid_hole_y = mean(base_y_coords);
    
    % Find nearest grid indices to centroid
    [~, center_x] = min(abs(x - centroid_hole_x));
    [~, center_y] = min(abs(y - centroid_hole_y));
    
    % If you want to ensure a true grid-center fallback if centroid is NaN/zero:
    if isnan(centroid_hole_x) || isnan(centroid_hole_y) || isempty(center_x) || isempty(center_y)
        center_x = round(Nx/2);
        center_y = round(Ny/2);
        warning('Horn base centroid invalid — falling back to grid center.');
    end
    
    
    % Create spatial grid for membrane
    [X_2d, Y_2d] = meshgrid(x, y);
    X_2d_centered = X_2d - x(center_x);
    Y_2d_centered = Y_2d - y(center_y);
    R_2d = sqrt(X_2d_centered.^2 + Y_2d_centered.^2);
    
    % Membrane mask (circular region)
    pmut_mask_2d = R_2d <= pmut_radius;
    
    % Velocity profile using Bessel function J0
    velocity_profile = zeros(Ny, Nx);
    velocity_profile(pmut_mask_2d) = v0 * besselj(0, k_membrane * R_2d(pmut_mask_2d));
    
    % Normalize to ensure peak velocity = v0
    velocity_profile = velocity_profile / max(abs(velocity_profile(:))) * v0;
    
    % Define the source as a velocity boundary condition
    source = struct();
    source.u_mask = zeros(Ny, Nx, Nz);
    source.u_mask(:,:,source_z_index) = pmut_mask_2d;
    
    % Flatten spatial mask for k-Wave source definition
    n_source_points = sum(pmut_mask_2d(:));
    nt = length(kgrid.t_array);
    
    % Create time-varying velocity signal
    source.uz = zeros(n_source_points, nt);
    for it = 1:nt
        source.uz(:, it) = velocity_profile(pmut_mask_2d) * sin(2*pi*f0*kgrid.t_array(it));
    end
    
    %% ========== SENSOR: PRESSURE MEASUREMENT PLANE ==========
  
    % Find indices closest to observation heights
    z_up_target_abs = z_base + z_up_target;        
    [~, z_obs_idx] = min(abs(z - z_up_target_abs));

    % Sensor records entire plane at observation height
    sensor.mask = zeros(Ny, Nx, Nz);
    sensor.mask(:,:,z_obs_idx) = 1;
    
    % sensor.mask(:,:,z_obs_idx_down) = 1;
    
    % Also record xz plane (y = center) for visualization
    y_center_idx = round(Ny/2);
    for iz = 1:Nz
        sensor.mask(y_center_idx, :, iz) = 1;
    end
    
    % Record maximum pressure at each point
    sensor.record = {'p_max'};
        
    %% ========== MEDIUM SETUP ==========
    medium.sound_speed = sound_speed;  % Uniform air
    medium.density = density;          % Uniform air
    
    % Horn region
    medium.sound_speed(horn_mask == 1) = c_horn;  % Keep same to avoid CFL issues
    medium.density(horn_mask == 1) = rho_horn;  % Much higher density = rigid wall
    
    %% ========== RUN k-WAVE SIMULATION ==========

    input_args = {'PMLSize', pml_size, 'DataCast', 'gpuArray-single', ...
                      'PlotSim', false, 'PlotPML', false};
    sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});

    %% ========== PROCESS RESULTS ==========
    
    % Extract pressure data
    p_max_data = gather(sensor_data.p_max);

    % Reconstruct pressure field from sensor data
    p_max_field = zeros(Ny, Nx, Nz);
    sensor_indices = find(sensor.mask);
    for i = 1:length(sensor_indices)
        p_max_field(sensor_indices(i)) = p_max_data(i);
    end
    
    % Extract observation plane (xy at z_obs_idx)
    p_max_2d = squeeze(p_max_field(:, :, z_obs_idx));
    
    % Extract xz plane (at y = center)
    p_max_xz = squeeze(p_max_field(y_center_idx, :, :))';
    
    % Find maximum pressure on observation plane
    [max_pressure, ~] = max(p_max_2d(:));
 
    fprintf('\n\n Amplitude improvement ratio: %.6f\n\n\n', ...
            (max_pressure)/(base_amplitude));
     
    improvement_ratio = (max_pressure)/(base_amplitude);

    % Apply penalty if improvement ratio is unrealistically large (unstable
    % model)
    if improvement_ratio < 5000
        cost = - (max_pressure)/(base_amplitude);
    else
        cost = 10000000;
    end
    
end
