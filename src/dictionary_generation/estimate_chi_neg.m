function chi_neg = estimate_chi_neg(circData, theta, g_ratio, chi_iso, chi_ani)
% ESTIMATE_CHI_NEG Estimates the negative susceptibility component.
% 
% This function calculates the expected negative susceptibility (chi_neg) 
% by modeling the myelin field perturbation and fitting it against a 
% theoretical point dipole field.
%
% Inputs:
%   circData : Matrix containing circle distributions (radii and centers)
%   theta    : Fiber orientation angle
%   g_ratio  : Ratio of inner to outer radius of the myelin sheath
%   chi_iso  : Isotropic susceptibility of myelin
%   chi_ani  : Anisotropic susceptibility of myelin
%
% Outputs:
%   chi_neg  : Estimated negative susceptibility component

%% 1. Initialize 3D Grid
matrix_size = 257;
lim = (matrix_size - 1) / 2;
x = -lim:lim;
y = -lim:lim;
z = -lim:lim;
[X,Y,Z] = meshgrid(x,y,z);

%% 2. Compute Point Dipole Field
dipole_location = [0 0 0];
dipole_susceptibility = 1; 
dipole_field = field_given_point_dipole(dipole_location,dipole_susceptibility,X,Y,Z);

% dipole_field((matrix_size+1)/2-2:(matrix_size+1)/2+2,(matrix_size+1)/2-2:(matrix_size+1)/2+2,(matrix_size+1)/2-2:(matrix_size+1)/2+2) = 0;
dipole_field(isnan(dipole_field)) = 0;

%% 3. Rotate and Extract Central Slice
direction = [1 0 0];
dipole_field = imrotate3(dipole_field, theta, direction,'linear','crop','FillValues',0);
dipole_field = squeeze(dipole_field(:,:,(matrix_size+1)/2));

%% 4. Calculate Myelin Field Perturbation
[phi,r] = polar_mesh(x,y);
radii = circData(:,3) / 512;
[counts, radii] = groupcounts(radii);

myelin_field = zeros(size(phi));

for i = 1:length(radii)
    r_o = radii(i);
    r_i = r_o * g_ratio;
    
    % Superimpose the field contribution from each hollow cylinder size
    myelin_field = myelin_field + counts(i) * ((chi_iso+chi_ani) .* sind(theta)^2 .* cosd(2*phi)) ./ 2 .* ((r_o^2 - r_i^2)./(r).^2);
    
    % Mask the origin to avoid singularity issues
    myelin_field((matrix_size+1)/2,(matrix_size+1)/2) = 0;
    % myelin_field((matrix_size+1)/2-2:(matrix_size+1)/2+2,(matrix_size+1)/2-2:(matrix_size+1)/2+2) = 0;
end
% end

%% 5. Fit Linear Model to Estimate chi_neg
fit = fitlm(dipole_field(:),myelin_field(:));
chi_neg = table2array(fit.Coefficients(2,1));

end

%% Local Helper Function
function field = field_given_point_dipole(dipole_location,dipole_susceptibility,X,Y,Z)
% Calculates the magnetic field perturbation of a point dipole
    field = dipole_susceptibility*(3*(Z-dipole_location(3)).^2 - ...
        (X-dipole_location(1)).^2-(Y-dipole_location(2)).^2-(Z-dipole_location(3)).^2) ...
        ./(4*pi*((X-dipole_location(1)).^2+(Y-dipole_location(2)).^2+(Z-dipole_location(3)).^2.^2).^2.5);
end


