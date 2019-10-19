clear; clc; close all;
LOAD_RESULTS = 0;

if LOAD_RESULTS
    %% Load pre-computed NR results
    load('results.mat');
    addpath(genpath('ui_utils\'));
    nbands = size(msi, 3);
else
    %% Run each method to obtain the NR results
    Comparison;
end

%% Diaplay QA indices
methods = {'Noisy', 'BwKSVD', 'BwBM3D', 'IntKSVD', '3DNLM', 'BM4D', 'LRTA', 'PARAFAC', 'TensorDL'};
fprintf('%10s%10s%10s%10s%10s%10s%10s\n', 'method', 'PSNR', 'SSIM', 'FSIM', 'ERGAS', 'SAM', 'time(s)');
fprintf('----------------------------------------------------------------------\n');
for i = 1:9
    if 0 == ENABLE_BITS(i)
        fprintf('%10s%10s%10s%10s%10s%10s%10s\n', ...
            methods{i}, '--', '--', '--', '--', '--', '--');
    else
        fprintf('%10s%10.2f%10.4f%10.4f%10.2f%10.4f%10.2f\n', ...
            methods{i}, psnr(i), ssim(i), fsim(i), ergas(i), sam(i), time(i));
    end
end

%% Show the MSIs
initUI;