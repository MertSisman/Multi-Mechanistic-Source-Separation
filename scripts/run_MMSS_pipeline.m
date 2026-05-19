clear
close all
clc
rng("default")

%% ========================================================================
%  USER CONFIGURATION
%  ========================================================================
%  Define the absolute path to the main MMSS repository folder.
%    (e.g., '/home/username/Documents/GitHub/MMSS')
REPO_PATH = '/path/to/MMSS_Repo';

% Define the path to your Python executable. 
% If using a specific Conda/Venv environment, provide the absolute path.
% (e.g., '/home/user/miniconda3/envs/MMSS_env/bin/python3')
% If Python is in your system PATH, you can simply use 'python3'.
PYTHON_EXEC = 'python3';
% =========================================================================
%  AUTOMATED PATH SETUP (Do not change)
% =========================================================================
% Define subdirectories based on standard repo architecture
DATA_DIR   = fullfile(REPO_PATH, 'data');
SCRIPT_DIR = fullfile(REPO_PATH, 'scripts');
DICTI_DIR   = fullfile(REPO_PATH, 'dictionary');
ML_DIR   = fullfile(REPO_PATH, 'src','ml');

PYTHON_SCRIPT_PATH = fullfile(ML_DIR, 'train_and_process.py');

%% 1. Load Simulated Dictionary and in vivo MRI Data
disp('Loading dictionary and in vivo data...');

% Note: Assuming the generated dictionary is stored in the data/ folder
load(fullfile(DICTI_DIR, 'dictionary_MISSS_100K.mat'))

% Optional: If you had a 'dirs' loop in a wrapper, ensure 'm' is defined. 
% disp(strcat("Starting ", dirs(m).name));

load(fullfile(DATA_DIR, 'iField.mat'), 'TE', 'iField')
load(fullfile(DATA_DIR, 'QSM.mat'))
load(fullfile(DATA_DIR, 'RDF.mat'), 'Mask')

%% 2. Dictionary Preparation for ML Training
L = length(dictionary.Afw);

dictionary.magnitude = abs(dictionary.signal);

TE_original =  dictionary.TE * 1e3;
TE_target = TE * 1e3;
N_TE = length(TE);

SNR0 = rand(L,1)*40+10;
SNR = repmat(SNR0,[1,N_TE]);

degree = 5;

dictionary.magnitude_interpolated = interpolate_dictionary(dictionary.magnitude, TE_original, TE_target, degree);
dictionary.magnitude_interpolated = dictionary.magnitude_interpolated + randn(size(dictionary.magnitude_interpolated))*1./SNR;
dictionary.magnitude_interpolated_normalized = dictionary.magnitude_interpolated./vecnorm(dictionary.magnitude_interpolated,2,2);

chi_ani = dictionary.chi_ani * 1e6; %ppm
chi_iron = dictionary.chi_iron * 1e6; %ppm
chi_iso = dictionary.chi_iso*1e6;  %ppm
MVF = dictionary.FVF.*(1-(dictionary.g_ratio).^2);
chi_neg = dictionary.chi_neg*1e6;
chi_pos = dictionary.chi_pos*1e6;

chi_total = chi_neg + chi_pos;
chi_total = chi_total + randn(size(chi_total))*0.3/SNR0.';
chi_total = (chi_total + 0.1) / 0.4;

magnitude = dictionary.magnitude_interpolated_normalized;

chi_neg = chi_neg*(-1) / 0.1;
chi_pos = chi_pos / 0.3;

theta = dictionary.theta;
theta = theta / 90;

fwf = dictionary.fwf;

data = [magnitude chi_total.'];
labels = [chi_neg; chi_pos; fwf; MVF; theta].';

% Save training data directly into the DATA_DIR
save(fullfile(DATA_DIR, 'training_data_MLP_with_QSM_modeling.mat'), 'data', 'labels', 'N_TE')

%% 3. In vivo Data Preparation for Inference
magnitude = abs(iField) .* Mask;
s = size(magnitude);
magnitude = magnitude ./ vecnorm(magnitude,2,4);
magnitude(isnan(magnitude)) = 0;
magnitude(isinf(magnitude)) = 0;

QSM = (QSM + 0.1) / 0.4;

input = cat(4,magnitude,QSM);
input = reshape(input,prod(size(input,1:3)),size(input,4));

% Save inference input directly into the DATA_DIR
save(fullfile(DATA_DIR, 'data.mat'), 'input', '-v7.3')

%% 4. Execute MLP via Python
disp('Executing PyTorch Multi Layer Perceptron...');

if ~isfile(PYTHON_SCRIPT_PATH)
    error('Cannot find the Python script at: %s', PYTHON_SCRIPT_PATH);
end

% Temporarily change to the DATA directory so the Python script 
% can load/save the .mat files natively without pathing modifications.
original_dir = pwd;
cd(DATA_DIR);

% Execute Python
command = sprintf('%s "%s"', PYTHON_EXEC, PYTHON_SCRIPT_PATH);
[status, cmdout] = system(command);

% Return to the original directory immediately after execution
cd(original_dir);

if status ~= 0
    error('Python execution failed with the following error:\n%s', cmdout);
else
    disp('Python inference completed successfully.');
end

%% 5. Reconstruct Output Maps
% Load the results from the DATA_DIR
load(fullfile(DATA_DIR, 'MMSS_results.mat'))

chi_neg = reshape(MMSS_results(:,1),size(QSM)) * 100 .*Mask;
chi_pos = reshape(MMSS_results(:,2),size(QSM)) * 300 .*Mask;
fwf = reshape(MMSS_results(:,3),size(QSM)) * 100 .*Mask;
MVF = reshape(MMSS_results(:,4),size(QSM)) * 100 .*Mask;
theta = reshape(MMSS_results(:,5),size(QSM)) * 90 .*Mask;

% Save final outputs back to DATA_DIR
save(fullfile(DATA_DIR, 'MMSS_results.mat'), 'chi_neg', 'chi_pos', 'fwf', 'MVF', 'theta')
disp('MMSS Pipeline complete. Results saved to data directory.');

%% 6. Cleanup Intermediate Files
disp('Cleaning up intermediate files...');

training_file = fullfile(DATA_DIR, 'training_data_MLP_with_QSM_modeling.mat');
if isfile(training_file)
    delete(training_file);
end

inference_file = fullfile(DATA_DIR, 'data.mat');
if isfile(inference_file)
    delete(inference_file);
end

disp('Cleanup complete.');


%% Local Helper Functions
function D_interpolated = interpolate_dictionary(D, TE_original, TE_target, degree)

if size(TE_original,2) ~= 1
    TE_original = TE_original';
end
if size(TE_target,2) ~= 1
    TE_target = TE_target';
end

A_original = [];
for i = 1:degree
    A_original = cat(2,A_original,TE_original.^(i-1));
end

A_target = [];
for i = 1:degree
    A_target = cat(2,A_target,TE_target.^(i-1));
end

D_log = log(D);

D_interpolated_log = (A_target * pinv(A_original) * D_log.').';

D_interpolated = exp(D_interpolated_log);
end