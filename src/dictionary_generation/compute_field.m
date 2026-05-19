function [field] = compute_field(matrix_size, r_o, g_ratio, theta, chi_iso, chi_ani)
% COMPUTE_FIELD Computes 2D field distribution given a hollow circular distribution.
%
% Derived in Wharton, S. and Bowtell, R., 2012. Fiber orientation-dependent 
% white matter contrast in gradient echo MRI. Proceedings of the National 
% Academy of Sciences, 109(45), pp.18559-18564.
%
% Inputs: 
%   matrix_size = Size of the square grid
%   r_o         = Outer radius including myelin sheath
%   g_ratio     = Ratio of inner radius to outer radius (r_i / r_o)
%   theta       = Fiber orientation angle relative to B0
%   chi_iso     = Isotropic susceptibility of myelin
%   chi_ani     = Anisotropic susceptibility of myelin
%
% Output: 
%   field       = Total 2D magnetic field perturbation

%% 1. Define Geometry and Coordinate System
r_i = ceil(r_o*g_ratio); % Inner radius excluding myelin.

% Generate Cartesian grid.
lim = (matrix_size - 1) / 2;
x = -lim:lim;
y = -lim:lim;

% Convert Cartesian grid to polar grid.
[phi,r] = polar_mesh(x,y);

% Create masks that will represent inner (intracellular), annular (myelin),
% and outer (extracellular) regions.
inner_region = double(abs(r) < r_i);
annular_region = double((abs(r) >= r_i) .* (abs(r) < r_o));
outer_region = double(abs(r) >= r_o);

% Rotate azimuthal angle to align correctly with the B0 orientation conventions
phi = imrotate(phi,90);

%% 2. Field due to Isotropic Susceptibility
% Calculate the isotropic perturbation in the 3 regions.
delta_f_outer_iso = (chi_iso .* sind(theta)^2 .* cosd(2*phi)) ./ 2 .* ((r_o^2 - r_i^2)./(r).^2);
delta_f_annular_iso = chi_iso ./ 2 .* (cosd(theta)^2 - 1/3 - sind(theta)^2 .* cosd(2*phi).* (r_i^2./r.^2));
delta_f_inner_iso = 0;

% Combine regions using spatial masks
delta_f_iso = delta_f_inner_iso .* inner_region + delta_f_annular_iso .* annular_region + delta_f_outer_iso .* outer_region;

%% 3. Field due to Anisotropic Susceptibility
% Calculate the anisotropic perturbation in the 3 regions.
delta_f_outer_ani = (chi_ani .* (sind(theta)).^2 .* cosd(2*phi))./8 .* ((r_o^2 - r_i^2)./r.^2);
delta_f_annular_ani = chi_ani .* (sind(theta).^2.*(-5/12 - cosd(2*phi)/8 .* (1+r_i^2./r.^2) + 3/4 .* log(r_o./r))-cosd(theta)^2/6);
delta_f_inner_ani = (3 .* chi_ani  .* sind(theta).^2) / 4 .* log(r_o./r_i);

% Combine regions using spatial masks
delta_f_ani = delta_f_inner_ani .* inner_region + delta_f_annular_ani .* annular_region + delta_f_outer_ani .* outer_region;

% Handle singularities/NaNs at the origin
delta_f_ani(isnan(delta_f_ani)) = 0;

%% 4. Total Field Perturbation
% Sum the isotropic and anisotropic components
field = delta_f_ani + delta_f_iso;

end