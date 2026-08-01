% Test code to validate, visualize, and extract dxf and STL file designs 
% produced by optimize_acoustic_housing.m
% Code by Mostafa Sedky
% August 1, 2026

% You can use this to establish a baseline acoustic pressure without a horn
% so you can plug it into acoustic_cost.m

clear; clc; close all;

%% ========== PARAMETERS ==========
% PMUT parameters
pmut_diameter = 0.7e-3;           % pMUT diameter (m)
pmut_radius = pmut_diameter/2;
f0 = 180e3;                     % Resonance frequency (Hz)
v0 = 1;                      % Peak particle velocity (m/s)

% Acoustic housing (horn) parameters
cylinder_diameter = 8e-3;       % Outer cylinder diameter (m)
cylinder_radius = cylinder_diameter/2;
cylinder_height = 0.475e-3;      % Base packaging cylinder height (m)
bottom_cylinder_radius = pmut_radius * 1.2;
makeHorn = 1;
    
% Generalized ellipsoidal hole shape parameters
horn_shape.r_base = pmut_radius * 2.4;  % Base radius at pMUT

% HORN OPTIMIZATION PARAMETERS:
horn_height = 0.0031;             % Height of acoustic housing (m)
horn_shape.r_top = 0.0018;                % Top radius
horn_shape.p = 4.9299;                        % Ellipsoid exponent p
horn_shape.q = 6.5188;                        % Ellipsoid exponent q
horn_shape.r = 2.3663;                        % Ellipsoid exponent r
horn_shape.A_theta = 0;       % Azimuthal modulation amplitude
horn_shape.m_theta = 0;         % Azimuthal mode number
horn_shape.A_x = 0;           % X-direction modulation
horn_shape.A_y = 0;           % Y-direction modulation
horn_shape.A_z = 0;           % Z-direction modulation
horn_shape.m_z = 0;             % Z-direction mode number

% Medium properties
c_air = 343;                    % Speed of sound in air (m/s)
rho_air = 1.21;                 % Density of air (kg/m^3)
c_horn = 2730;                  % Speed of sound in horn material (m/s) - e.g., aluminum
rho_horn = 1190;                % Density of horn material (kg/m^3)

% Simulation parameters
z_up_target = 4e-3;             % Target height to measure pressue (m)
z_down_target = 4e-3;           % Target height to measure pressure (m)(use if one needs to check the pressure at two different heights)
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

fprintf('Grid size: %d x %d x %d = %.2f M points\n', Nx, Ny, Nz, Nx*Ny*Nz/1e6);
fprintf('Spatial resolution: dx = %.4f mm (lambda/%.1f)\n', dx*1e3, lambda/dx);

kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

% Time array
kgrid.Nt = Nsteps;
kgrid.dt = deltaT;

fprintf('Time steps: %d\n', length(kgrid.t_array));
fprintf('dt = %.2e s (CFL ~ %.2f)\n', kgrid.dt, c_air*kgrid.dt/dx);
fprintf('\n\n');

%% ========== GENERATE ACOUSTIC HOUSING GEOMETRY ==========
fprintf('Generating acoustic housing geometry...\n');

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
scale = @(idx) (z_horn(idx) / horn_height); %.^horn_shape.r;  % Normalized height
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

fprintf('Horn boundary mask created: %d solid points\n', sum(horn_mask(:)));
fprintf('Acoustic housing generated.\n');
fprintf('\n\n');

%% ===== EXPORT SMOOTH STL SURFACE OF ACOUSTIC HORN IN MM =====

fprintf('Exporting STL surface in millimeters...\n');

% Scale factor (meters to millimeters):
scaleF = 1000;

% Shift Z to correct physical location and scale all coordinates
Zs = (horn_Z + z_base) * scaleF;
Xs = horn_X * scaleF;
Ys = horn_Y * scaleF;

% Convert parametric surface into faces + vertices
[F, V] = surf2patch(Xs, Ys, Zs, 'triangles');

% Build triangulation object (required by your stlwrite)
TR = triangulation(F, V);

% Write STL file in mm
stlwrite(TR, 'figure_1_acoustic_horn_mm.stl');

fprintf('STL file exported: acoustic_horn_mm.stl (units in mm)\n');

%% ===== EXPORT DXF OF AXISYMMETRIC PROFILE WITH ORIGIN AT MIN-Z =====

% Raw profile in mm
z_profile = (z_base + z_horn) * 1000;   % mm
r_profile = r_curve * 1000;             % mm

% Shift so minimum z is the origin
z0 = min(z_profile);
z_shifted = z_profile - z0;

filename = 'figure_1_horn_profile_origin_at_base.dxf';
fid = fopen(filename, 'w');

fprintf(fid, '0\nSECTION\n2\nENTITIES\n');

for i = 1:length(z_shifted)-1
    fprintf(fid, '0\nLINE\n8\n0\n');
    fprintf(fid, '10\n%f\n20\n%f\n30\n0\n', z_shifted(i), r_profile(i));
    fprintf(fid, '11\n%f\n21\n%f\n31\n0\n', z_shifted(i+1), r_profile(i+1));
end

fprintf(fid, '0\nENDSEC\n0\nEOF\n');
fclose(fid);

fprintf('DXF curve exported with base at (0,0): %s\n', filename);

%% ========== SOURCE: pMUT VELOCITY PROFILE ==========
fprintf('Creating pMUT velocity source...\n');

% For a circular clamped membrane vibrating in (0,1) mode:
% v(r,t) = v0 * J0(k*r) * sin(omega*t)

% Membrane parameters
k_membrane = 2.4048 / pmut_radius;  % Wavenumber for (0,1) mode

% --- Find horn base centroid (use horn profile at z = 0) ---
% NOTE: horn_X, horn_Y were created as arrays indexed by z_horn index.
% Use the first z_horn row which corresponds to base (z=0).
base_x_coords = horn_X(1,:);   % coordinates of hole boundary in x (m)
base_y_coords = horn_Y(1,:);   % coordinates of hole boundary in y (m)

% Compute centroid of the hole boundary (in physical coords)
centroid_hole_x = mean(base_x_coords);
centroid_hole_y = mean(base_y_coords);

% Find nearest grid indices to centroid
[~, center_x] = min(abs(x - centroid_hole_x));
[~, center_y] = min(abs(y - centroid_hole_y));

% True grid-center fallback if centroid is NaN/zero:
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
fprintf('Creating time-varying velocity signal (%d points, %d timesteps)...\n', ...
        n_source_points, nt);

source.uz = zeros(n_source_points, nt);
for it = 1:nt
    source.uz(:, it) = velocity_profile(pmut_mask_2d) * sin(2*pi*f0*kgrid.t_array(it));
end

fprintf('pMUT membrane velocity source created at z = %.2f mm\n', z(source_z_index)*1e3);
fprintf('Membrane mode: (0,1) with Bessel J0 profile\n');
fprintf('Peak velocity: %.3f mm/s\n', v0*1e3);
fprintf('Velocity at center: %.3f mm/s\n', velocity_profile(center_y, center_x)*1e3);
fprintf('Velocity at edge: %.3f mm/s\n', velocity_profile(round(center_y + pmut_radius/dy), center_x)*1e3);
fprintf('\n\n');

%% ========== SENSOR: PRESSURE MEASUREMENT PLANE ==========
fprintf('Setting up sensor plane...\n');

% Find indices closest to observation heights
z_up_target_abs = z_base + z_up_target;        % if z_up_target is relative
z_down_target_abs = z_base + z_down_target;
[~, z_obs_idx] = min(abs(z - z_up_target_abs));
[~, z_obs_idx_down] = min(abs(z - z_down_target_abs));

fprintf('Top observation plane at z = %.2f mm (index %d)\n', z(z_obs_idx)*1e3, z_obs_idx);
fprintf('Bottom observation plane at z = %.2f mm (index %d)\n', z(z_obs_idx_down)*1e3, z_obs_idx_down);

% Sensor records entire plane at observation height
sensor.mask = zeros(Ny, Nx, Nz);
sensor.mask(:,:,z_obs_idx) = 1;

sensor.mask(:,:,z_obs_idx_down) = 1;

% Also record xz plane (y = center) for visualization
y_center_idx = round(Ny/2);
for iz = 1:Nz
    sensor.mask(y_center_idx, :, iz) = 1;
end

% Record maximum pressure at each point
sensor.record = {'p_max'};

fprintf('Recording: xy plane at z=%.1f mm and xz plane at y=0\n', z(z_obs_idx)*1e3);
fprintf('\n\n');

%% ========== MEDIUM SETUP ==========
medium.sound_speed = sound_speed;  % Uniform air
medium.density = density;          % Uniform air

% % Air absorbption properties (did not need, but can add for some cases)
% medium.alpha_coeff = 1.8;  % dB/(MHz^y cm)
% medium.alpha_power = 1.0;  % Exponent 

% Horn region
medium.sound_speed(horn_mask == 1) = c_horn; 
medium.density(horn_mask == 1) = rho_horn;  

% % Horn absorbption properties (did not need, but can add for some cases)
% medium.alpha_coeff(horn_mask == 1) = 0;  % dB/(MHz^y cm)

fprintf('\n=== Medium Properties ===\n');
fprintf('Air regions: c = %.1f m/s, rho = %.2f kg/m^3\n', c_air, rho_air);
fprintf('Horn boundary regions: c = %.1f m/s, rho = %.1f kg/m^3 (rigid)\n', ...
        c_horn, rho_horn);
fprintf('Grid spacing: dx = %.4f mm\n', dx*1e3);
fprintf('Frequency: f0 = %.1f kHz\n', f0/1e3);
fprintf('Wavelength in air: %.2f mm (%.1f ppw)\n', lambda*1e3, lambda/dx);


%% ========== RUN k-WAVE SIMULATION ==========
fprintf('\nStarting k-Wave simulation...\n');
tic;

% Use GPU if available (note: this only works for NVIDIA GPUs as it uses CUDA)
try
    gpuDevice(1);
    use_gpu = true;
    fprintf('GPU detected - using GPU acceleration\n');
    input_args = {'PMLSize', pml_size, 'DataCast', 'gpuArray-single', ...
                  'PlotSim', false, 'PlotPML', false};
catch
    use_gpu = false;
    fprintf('No GPU detected - using CPU\n');
    input_args = {'PMLSize', pml_size, 'PlotSim', false, 'PlotPML', false};
end

% Check source magnitude before running
fprintf('\n=== Source Check ===\n');
fprintf('Max source pressure: %.2f Pa\n', max(abs(source.uz(:))));
fprintf('Number of source points: %d\n', n_source_points);
fprintf('Source area: %.2f mm^2\n', n_source_points * dx * dy * 1e6);
fprintf('Estimated acoustic pressure amplitude: %.2f Pa (assuming rho*c*v)\n', ...
        rho_air * c_air * max(abs(source.uz(:))));

sensor_data = kspaceFirstOrder3DG(kgrid, medium, source, sensor, input_args{:});

t_sim = toc;
fprintf('Simulation completed in %.2f seconds\n', t_sim);
%% ========== PROCESS RESULTS ==========
fprintf('\nProcessing results...\n');

% Extract pressure data
if use_gpu
    p_max_data = gather(sensor_data.p_max);
else
    p_max_data = sensor_data.p_max;
end

% Reconstruct pressure field from sensor data
p_max_field = zeros(Ny, Nx, Nz);
sensor_indices = find(sensor.mask);
for i = 1:length(sensor_indices)
    p_max_field(sensor_indices(i)) = p_max_data(i);
end

% Extract observation plane (xy at z_obs_idx)
p_max_2d = squeeze(p_max_field(:, :, z_obs_idx));

% Extract observation plane (xy at z_obs_idx_down)
p_max_2d_down = squeeze(p_max_field(:, :, z_obs_idx_down));

% Extract xz plane (at y = center)
p_max_xz = squeeze(p_max_field(y_center_idx, :, :))';

% Find maximum pressure on observation plane
[max_pressure, max_idx] = max(p_max_2d(:));
[max_y, max_x] = ind2sub([Ny, Nx], max_idx);
max_x_pos = x(max_x);
max_y_pos = y(max_y);

[max_pressure_down, max_idx_down] = max(p_max_2d_down(:));
[max_y_down, max_x_down] = ind2sub([Ny, Nx], max_idx_down);
max_x_pos_down = x(max_x_down);
max_y_pos_down = y(max_y_down);

fprintf('\n========== RESULTS ==========\n');
fprintf('Maximum pressure up: %.2f Pa (%.2f dB re 20 µPa)\n', ...
        max_pressure, 20*log10(max_pressure/20e-6));
fprintf('Maximum pressure down: %.2f uPa (%.2f dB re 20 µPa)\n', ...
        max_pressure_down * 1e6, 20*log10(max_pressure_down/20e-6));
fprintf('Amplitude improvement (top/bottom): %.2f\n', ...
        max_pressure/max_pressure_down);
fprintf('Location: x = %.2f mm, y = %.2f mm\n', max_x_pos*1e3, max_y_pos*1e3);
fprintf('At top height: z = %.2f mm and bottom height: z = %.2f mm\n', z(z_obs_idx)*1e3, z(z_obs_idx_down)*1e3);

% Calculate on-axis pressure
p_center = p_max_2d(round(Ny/2), round(Nx/2));
fprintf('On-axis pressure: %.2f Pa (%.2f dB re 20 µPa)\n', ...
        p_center, 20*log10(p_center/20e-6));

%% ========== REMOVE PML FROM SLICES BEFORE PLOTTING ==========
% Define crop indices that exclude PML layers at all sides
ix_crop = (pml_size+1):(Nx-pml_size);
iy_crop = (pml_size+1):(Ny-pml_size);
iz_crop = (pml_size+1):(Nz-pml_size);

% Cropped coordinate vectors
x_crop = x(ix_crop);
y_crop = y(iy_crop);
z_crop = z(iz_crop);

%% ========== VISUALIZATION ==========

%% Figure 1: Pressure distributions
figure(); %'Position', [100, 100, 1400, 500]

% Plot 1: Pressure distribution at observation plane 1 (xy)
subplot(1,4,1);
imagesc(x*1e3, y*1e3, p_max_2d(ix_crop, iy_crop));
axis equal tight;
xlabel('x (mm)');
ylabel('y (mm)');
title(sprintf('Pressure at z = %.1f cm', z_up_target*1e2));
colorbar;
colormap(hot);

% Plot 2: Pressure distribution at observation plane 2 (xy)
subplot(1,4,2);
imagesc(x*1e3, y*1e3, p_max_2d_down(ix_crop, iy_crop));
axis equal tight;
xlabel('x (mm)');
ylabel('y (mm)');
title(sprintf('Pressure at z = %.1f cm', z_down_target*1e2));
colorbar;
colormap(hot);

% Plot 3: Horn geometry (analytical profile)
subplot(1,4,3);
hold on;
for iz = 1:min(n_z_horn, Nz)
    if z_horn(iz) <= horn_height
        plot(horn_X(iz,:)*1e3, z_horn(iz)*ones(size(horn_X(iz,:)))*1e3, 'b.', 'MarkerSize', 2);
    end
end
plot([-cylinder_radius, cylinder_radius]*1e3, [horn_height, horn_height]*1e3, 'k-', 'LineWidth', 2);
xlabel('Radius (mm)');
ylabel('Height (mm)');
title('Horn Profile (analytical)');
grid on;
axis equal;
ylim([0, horn_height*1e3*1.1]);

% Plot 4: Pressure distribution at observation plane (xz)
subplot(1,4,4);
imagesc(x*1e3, z*1e3, p_max_xz(iz_crop, ix_crop));
set(gca, 'YDir', 'normal');
axis equal tight;
xlabel('x (mm)');
ylabel('z (mm)');
title(sprintf('Pressure at xz plane'));
colorbar;
colormap(jet);


%% Figure 2: Voxelized horn domain and pressure in xz plane
figure(); % , [100, 650, 1400, 600]

% Plot 1: Voxelized horn domain (xz plane, y = 0)
y_center_idx = round(Ny/2);
xz_slice_material = squeeze(medium.sound_speed(:, y_center_idx, :))';  % Transpose for correct orientation
imagesc(x*1e3, z*1e3, xz_slice_material);
axis xy;
xlabel('x (mm)');
ylabel('z (mm)');
title('Voxelized Horn Domain (xz plane, y=0)');
colorbar;
colormap(gca, 'gray');
caxis auto;
hold on;
% Mark pMUT location
plot([x(center_x)-pmut_radius, x(center_x)+pmut_radius]*1e3, [z(source_z_index), z(source_z_index)]*1e3, 'r-', 'LineWidth', 3);
% Mark observation plane
plot([x(1), x(end)]*1e3, [z(z_obs_idx), z(z_obs_idx)]*1e3, 'g--', 'LineWidth', 2);
legend('pMUT', 'Observation plane', 'Location', 'best');


%% Figure 3: pMUT velocity profile
figure(); %, [100, 1300, 700, 600]

subplot(2,1,1);
imagesc(x*1e3, y*1e3, velocity_profile);
axis equal tight;
xlabel('x (mm)');
ylabel('y (mm)');
title('pMUT Membrane Velocity Profile (Bessel J_0)');
colorbar;
colormap(jet);
hold on;
% Draw membrane boundary
theta_circle = linspace(0, 2*pi, 100);
plot(x(center_x)*1e3 + pmut_radius*1e3*cos(theta_circle), ...
     y(center_y)*1e3 + pmut_radius*1e3*sin(theta_circle), 'w--', 'LineWidth', 2);

subplot(2,1,2);
r_profile = linspace(0, pmut_radius, 100);
v_profile = v0 * besselj(0, k_membrane * r_profile);
v_profile = v_profile / max(abs(v_profile(:))) * v0;
plot(r_profile*1e3, v_profile*1e3, 'b-', 'LineWidth', 2);
xlabel('Radius (mm)');
ylabel('Velocity (mm/s)');
title('Radial Velocity Profile');
grid on;
hold on;
plot([0, pmut_radius]*1e3, [v0, v0]*1e3, 'r--', 'LineWidth', 1.5);
legend('J_0 profile', sprintf('Peak (%.2f mm/s)', v0*1e3), 'Location', 'best');


%% Figure 4: Enhanced xz pressure distribution with smooth horn overlay
figure();

% Plot pressure distribution using full coordinates
imagesc(x*1e3, z*1e3, p_max_xz);
set(gca, 'YDir', 'normal');
axis equal tight;
xlabel('x (mm)', 'FontSize', 12);
ylabel('z (mm)', 'FontSize', 12);
title('Pressure Distribution with Horn Geometry (xz plane, y=0)', 'FontSize', 14);
cb = colorbar;
ylabel(cb, 'Pressure (Pa)', 'FontSize', 11);
colormap(jet);
hold on;

% Overlay smooth horn profile boundary (right side)
horn_z_plot = (z_base + z_horn) * 1e3;  % Convert to mm
horn_r_plot = r_curve * 1e3;            % Convert to mm
plot(horn_r_plot, horn_z_plot, 'w-', 'LineWidth', 2.5);

% Overlay smooth horn profile boundary (left side - mirror)
plot(-horn_r_plot, horn_z_plot, 'w-', 'LineWidth', 2.5);

% Draw outer cylinder walls
plot([cylinder_radius, cylinder_radius]*1e3, [z_base*1e3, (z_base + horn_height)*1e3], 'w-', 'LineWidth', 2.5);
plot([-cylinder_radius, -cylinder_radius]*1e3, [z_base*1e3, (z_base + horn_height)*1e3], 'w-', 'LineWidth', 2.5);

% Draw bottom cylinder below horn
z_bottom = z_base * 1e3;
z_cyl_bottom = (z_base - cylinder_height) * 1e3;
plot([bottom_cylinder_radius, bottom_cylinder_radius]*1e3, [z_cyl_bottom, z_bottom], 'w-', 'LineWidth', 2.5);
plot([-bottom_cylinder_radius, -bottom_cylinder_radius]*1e3, [z_cyl_bottom, z_bottom], 'w-', 'LineWidth', 2.5);

% Draw horizontal connector if bottom cylinder is smaller than horn base
if bottom_cylinder_radius < horn_shape.r_base
    plot([bottom_cylinder_radius*1e3, horn_shape.r_base*1e3], [z_bottom, z_bottom], 'w-', 'LineWidth', 2.5);
    plot([-bottom_cylinder_radius*1e3, -horn_shape.r_base*1e3], [z_bottom, z_bottom], 'w-', 'LineWidth', 2.5);
end

% Mark observation plane (top)
plot([x(1), x(end)]*1e3, [z(z_obs_idx), z(z_obs_idx)]*1e3, 'g--', 'LineWidth', 2);

% Set x-axis limits
xlim([-2.7, 2.6]);
ylim([-2.4, 2.3])

fprintf('\nVisualization complete.\n');