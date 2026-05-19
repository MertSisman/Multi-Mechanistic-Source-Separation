clear 
close all
clc

%% Dictionary Generation for Multi Mechanistic Source Separation (MMSS)
% Mert Sisman 8/26/2024
% This script simulates a dictionary of gradient echo signals based on 
% microstructure informed susceptibility models.

rng('default') % for reproducibility 
matrix_size = 512;

%% 1. Set Volume Susceptibilities
chi_iso = -0.1e-6;    % -0.1 ppm / isotropic susceptibility of myelin sheaths.
chi_ani = -0.1e-6;    % -0.1 ppm / anisotropic susceptibility of myelin sheaths.
chi_iron = 0.3e-6;    %  0.3 ppm / isotropic iron susceptibility.
m = 1;

%% 2. Generate 2D Circle Distributions (Fiber Bundles)
for density = linspace(1e-3,1,10) % 2DFD
    
    % Define input parameters for 2D circle distribution generation.
    S.circSize = [matrix_size/180 matrix_size/160 matrix_size/140 matrix_size/120 matrix_size/100 matrix_size/80 matrix_size/60 matrix_size/40 matrix_size/24 matrix_size/20]; 
    S.nSizes = NaN;
    S.frameSize = [matrix_size*3/4 matrix_size*3/4];
    S.edgeType = 1;
    S.supressWarning = true;
    S.density = density;
    S.drawFrame = false;
    
    % Circle distribution generation.
    [circData, circHandles, frame, S] = bubblebath(S);
    circData = sortrows(circData,3);
    
    % Create a 2D grid.
    lim = (matrix_size - 1) / 2;
    x = -lim:lim;
    y = -lim:lim;
    [X,Y] = meshgrid(x,y);
    Y = flipud(Y);
    
    % Ensure that circle distribution (representing a fiber bundle)
    % resides entirely inside a circular region.
    for l = length(circData):-1:1
        if circData(l,1)^2 + circData(l,2)^2 > (matrix_size*3/8+circData(l,3))^2 
            circData(l,:) = []; 
        end    
    end
    
    %% 3. Monte Carlo Parameter Selection
    for dummy1 = 1:1000 % Randomly choose simulation input parameters 1000 times for each 2DFD.
        g_ratio = rand * 0.5 + 0.5;  % g-ratio (ri/ro).
        theta =  rand * 90;          % Fiber orientation.
        EID = rand;                  % Extracellular iron density.
        
        % Estimate chi negative.
        chi_neg = estimate_chi_neg(circData, theta, g_ratio, chi_iso, chi_ani);
        
        % Initialize masks for fibers (circles) and myelin (hollow circles) in 2D.
        circles = zeros(matrix_size);
        hollow_circles = zeros(matrix_size);
        
        % Initialize circle center mask.
        centers = zeros(matrix_size,matrix_size,length(S.circSize));
        
        for i = 1:length(circData)
            % Create a mask for the inner and outer circles of a single element.
            one_circle_outer = ((X-circData(i,1)).^2 + (Y-circData(i,2)).^2) < circData(i,3).^2;
            one_circle_inner = ((X-circData(i,1)).^2 + (Y-circData(i,2)).^2) < (circData(i,3).*g_ratio).^2;
            
            % Add the new circle to the total mask.
            circles = circles + one_circle_outer;
            hollow_circles = hollow_circles + one_circle_outer - one_circle_inner;
            
            % Determine the center of each circle and place it in the discrete center mask.
            for j = 1:length(S.circSize)  
                if circData(i,3) == S.circSize(j)
                    centers(round(circData(i,1)+matrix_size/2),round(circData(i,2)+matrix_size/2),j) = 1;
                end
            end
        end
        
        %% 4. Compute Field Perturbations
        % Initialize 2D myelin related field distribution.
        fields_multi_circle = zeros(matrix_size,matrix_size,length(S.circSize));
        
        % Compute the field perturbation due to each hollow circle and shift it
        % according to the center location.
        for j = 1:length(S.circSize)  
            field = compute_field(matrix_size,S.circSize(j),g_ratio, theta, chi_iso, chi_ani);
            fields_multi_circle(:,:,j) = conv2(centers(:,:,j),field,'same');
        end
        
        % Sum the field perturbations of each individual hollow circle.
        total_field_2D = sum(fields_multi_circle,3);
        
        %% 5. Extension to 3D and Iron Incorporation
        % Replicate 2D field, hollow circle mask, and circle mask to 3D.
        total_field_3D = repmat(total_field_2D, 1, 1, matrix_size); 
        myelin_mask = repmat(hollow_circles, 1, 1, matrix_size); 
        fiber_mask = repmat(circles, 1, 1, matrix_size); 
        
        % Resample the 3D field and the masks to the defined fiber orientation around -y axis.
        direction = [0 -1 0];
        total_field_3D = imrotate3(total_field_3D, theta, direction,'linear','crop','FillValues',0);
        myelin_mask = imrotate3(myelin_mask, theta, direction,'nearest','crop','FillValues',0);
        fiber_mask = imrotate3(fiber_mask, theta, direction,'nearest','crop','FillValues',0);
        
        % Extract the central ROI where the mGRE signal will be calculated. Crop the outer regions.
        myelin_field = total_field_3D(matrix_size/4+1:end-matrix_size/4,matrix_size/4+1:end-matrix_size/4,matrix_size/4+1:end-matrix_size/4);
        myelin_mask = myelin_mask(matrix_size/4+1:end-matrix_size/4,matrix_size/4+1:end-matrix_size/4,matrix_size/4+1:end-matrix_size/4);
        fiber_mask = fiber_mask(matrix_size/4+1:end-matrix_size/4,matrix_size/4+1:end-matrix_size/4,matrix_size/4+1:end-matrix_size/4);
        
        % Generate the extracellular iron distribution and the corresponding field.
        iron_volume = generate_iron_volume(fiber_mask,EID);
        iron_field = chi_iron * generate_iron_field(iron_volume);
        
        % Compute the total field.
        field = myelin_field + iron_field;
        
        % Calculate fiber volume fraction (FVF) and iron volume fraction (IVF).
        FVF = sum(fiber_mask(:)) / (matrix_size/2)^3;
        IVF = sum(iron_volume(:)) / (matrix_size/2)^3;
        
        %% 6. Calculate mGRE Signal
        % Define echo times.
        TE = 0:3e-3:60e-3;
        % Gyromagnetic ratio.
        gamma = 2 * pi * 42.577478518e6;
        % Scanner magnetic field.
        B_0 = 3;
        
        % T2 relaxation times
        T2_iew = 70e-3;     % Intra/extracellular water T2.
        T2_myelin = 16e-3;  % Myelin water T2.
        T2_FW= 2000e-3;      % Free Water T2.
        
        % Relative myelin proton density.
        relative_protondensity_myelin = 0.5;
        
        % Intra-extracellular water and myelin water weighted masks.
        iew_mask = 1 - myelin_mask;
        myelin_mask = relative_protondensity_myelin * myelin_mask;
        
        % Initialize mGRE signal.
        signal = zeros(1,length(TE));
        
        % Compute the T2 decays.
        T2decay_iew = exp(-TE/T2_iew);
        T2decay_myelin = exp(-TE/T2_myelin);
        T2decay_FW = exp(-TE/T2_FW);
        signal_iew = zeros(1,length(TE));
        signal_myelin = zeros(1,length(TE));
        
        % Compute complex factor distribution and sum it to get mGRE signal at
        % TE for myelin compartment and intra-extracellular compartment separately.
        for t = 1:length(TE)
            complex_factor_distribution_iew = iew_mask.*exp(-1i * TE(t) * gamma * B_0 * field );
            signal_iew(t) = sum(complex_factor_distribution_iew(:));
            
            complex_factor_distribution_myelin = myelin_mask.*exp(-1i * TE(t) * gamma * B_0 * field );
            signal_myelin(t) = sum(complex_factor_distribution_myelin(:));
        end
        
        % Calculate myelin and outside signal fractions.
        myelin_fraction = signal_myelin(1) / (signal_myelin(1) + signal_iew(1)); 
        iew_fraction = signal_iew(1) / (signal_myelin(1) + signal_iew(1)); 
        
        % Calculate signal excluding free water and normalize.
        signal_nonfw = T2decay_iew .* signal_iew + T2decay_myelin .* signal_myelin;
        signal_nonfw = signal_nonfw / abs(signal_nonfw(1));
        
        %% 7. Include Free Water Compartment
        % Dummy variable just a counter for the number of different free water fraction (fwf) values simulated.
        for dummy2 = 1:10
            fwf = rand;
            
            % Add the free water signal.
            signal = signal_nonfw * (1-fwf) + T2decay_FW * fwf;
            
            % Normalize the signal with the first echo.
            signal = signal / abs(signal(1));
            signal_norm = signal / norm(abs(signal));
            
            % Fit for R2*.
            f = fit(TE.',-log(abs(signal.')),'poly1');
            R2s = f.p1;
            
            % Calculate signal (water) fractions of different compartments.
            Amw = myelin_fraction*(1-fwf);
            Aiew = iew_fraction*(1-fwf);
            Afw = fwf;
            
            % Calculate free water volume fraction.
            FW_vf = Afw / (Afw + Aiew + Amw / relative_protondensity_myelin);
            
            %% 8. Append to Dictionary
            dictionary.R2s(m) = R2s;
            dictionary.FVF(m) = FVF*(1-FW_vf);
            dictionary.g_ratio(m) = g_ratio;
            dictionary.theta(m) = theta;
            dictionary.IVF(m) = IVF*(1-FW_vf);
            dictionary.fwf(m) = fwf;
            dictionary.chi_neg(m) = chi_neg*(1-FW_vf);
            dictionary.chi_pos(m) = IVF*chi_iron*(1-FW_vf);
            dictionary.Amw(m) = Amw;
            dictionary.Aiew(m) = Aiew;
            dictionary.Afw(m) = Afw;
            dictionary.signal(m,:) = signal;
            m = m+1;
        end
    end
    save dictionary_MMSS_100K.mat 'dictionary'
end

%% 9. Save Global Parameters
dictionary.chi_ani = chi_ani;
dictionary.chi_iso = chi_iso;
dictionary.chi_iron = chi_iron;
dictionary.TE = TE;
save dictionary_MMSS_100K.mat 'dictionary'

