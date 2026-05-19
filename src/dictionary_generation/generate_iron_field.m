function [iron_field] = generate_iron_field(iron_volume)
% GENERATE_IRON_FIELD Computes the magnetic field perturbation from iron.
%
% This function calculates the resulting magnetic field perturbation 
% produced by the generated iron volume. It uses a k-space approach, 
% multiplying the Fourier-transformed susceptibility by a dipole kernel, 
% and includes zero-padding to prevent circular convolution artifacts.
%
% Inputs:
%   iron_volume : 3D binary mask of iron locations
%
% Outputs:
%   iron_field  : 3D map of the magnetic field perturbation induced by the iron
%
% Dependencies:
%   dipole_kernel() - Helper function to generate the k-space dipole

%% 1. Get Matrix Dimensions
s = size(iron_volume);

%% 2. Generate Dipole Kernel
% Generate a dipole kernel twice the size of the volume for the field calculation
D = dipole_kernel(s*2, [1,1,1], [0,0,1]);

%% 3. Pad Susceptibility Matrix
% Zero-pad the susceptibility map to match the expanded kernel size
susceptibility_3D = padarray(iron_volume, s/2, 0, 'both');

%% 4. Calculate Field via FFT
% Perform k-space convolution (multiplication) and transform back to image space
iron_field = real(ifftn(fftn(susceptibility_3D) .* D));

%% 5. Crop to Original Size
% Crop the resulting field map back to the original input matrix dimensions
iron_field = iron_field(s(1)/2+1:end-s(1)/2, s(2)/2+1:end-s(2)/2, s(3)/2+1:end-s(3)/2);

end