function [iron_volume] = generate_iron_volume(fiber_volume, desired_iron_density)
% GENERATE_IRON_VOLUME Distributes extracellular iron randomly.
%
% This function generates a mask representing the spatial distribution 
% of iron particles. It ensures that iron is restricted entirely to the 
% extracellular space (i.e., outside the defined fiber volume).
%
% Inputs:
%   fiber_volume         : 3D mask of the fiber regions (1 = fiber, 0 = extracellular)
%   desired_iron_density : Probability density of iron occurrence (scalar between 0 and 1)
%
% Outputs:
%   iron_volume          : 3D binary mask of localized iron deposition

%% 1. Get Matrix Dimensions
s = size(fiber_volume);

%% 2. Distribute Iron
% Generate random iron locations strictly outside the fibers
iron_volume = (rand(s) < desired_iron_density) .* (1 - fiber_volume);

end