function [phi, r] = polar_mesh(x, y)
% POLAR_MESH Creates a 2D polar coordinate meshgrid.
%
% This helper function takes 1D Cartesian coordinate vectors and 
% generates 2D matrices for polar angle (in degrees) and radial distance.
% It ensures the Y-axis is flipped to match standard image coordinate systems.
%
% Inputs:
%   x : 1D array of x-coordinates
%   y : 1D array of y-coordinates
%
% Outputs:
%   phi : 2D matrix of polar angles (in degrees)
%   r   : 2D matrix of radial distances

%% 1. Create Cartesian Meshgrid
[X, Y] = meshgrid(x, y);

% Flip Y to match the expected image/matrix coordinate system
Y = flipud(Y);

%% 2. Convert to Polar Coordinates
% Note: cart2pol takes (Y, X) here based on the desired orientation
[phi, r] = cart2pol(Y, X);

%% 3. Convert Angle to Degrees
phi = rad2deg(phi);

end