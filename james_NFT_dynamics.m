clearvars;         
close all force;
addpath(genpath('functions_matlab'));

surface_interest = 'fsLR_32k';
hemisphere = 'lh';

% Load cortex mask
cortex = dlmread(sprintf('data/template_surfaces_volumes/%s_cortex-%s_mask.txt', surface_interest, hemisphere));

mesh_interest = 'midthickness';
[vertices, faces] = read_vtk(sprintf('data/template_surfaces_volumes/%s_%s-%s.vtk', surface_interest, mesh_interest, hemisphere));
surface_midthickness.vertices = vertices';
surface_midthickness.faces = faces';

num_modes = 200;

% =========================================================================
%         Simulation   
% =========================================================================

tmax = 864; % in s
tstep = 0.72;

param = loadParameters_wave_func;
param.tstep = tstep; % in s
param.tmax = tmax;  
param.tspan = [0, param.tmax];
param.T = 0:param.tstep:param.tmax;
param.is_time_ms = 0;

bal_param = loadParameters_balloon_func; % in s
bal_param.tstep = tstep;
bal_param.tmax = tmax;
bal_param.tspan = [0, param.tmax];
bal_param.T = 0:param.tstep:param.tmax;


% method = 'ODE';
method = 'Fourier';

param.r_s = 28;      % (default) in mm
param.gamma_s = 116; % (default) in s^-1
if param.is_time_ms==1
    param.gamma_s = 116*1e-3;
end

input = 'noise';
switch input
    case 'noise'
        seed = 1;
        rng(seed);
        mu = 0;
        std = 1;
        ext_input = randn(length(cortex), length(param.T)) * std + mu;
        str = sprintf("seed%d_mean%d_std%.2f", seed, mu, std);
    case 'pulse'
        cite = 15000;
        amp = 100;
        mu = 0;
        ext_input = ones(length(cortex), length(param.T)) * mu;
        ext_input(cite, 1:10) = ext_input(cite, 1:10) + amp;
        str = sprintf("cite%d_amp%d_mean%d", cite, amp, mu);
end

name = sprintf("james_%s_%s_%s_modes%d_T%d_tstep%.3f", input, str, method, num_modes, param.tmax, param.tstep);


% build the filename
outputDir = fullfile('.', 'james_mats');
if ~exist(outputDir, 'dir'), mkdir(outputDir); end
matFile = fullfile(outputDir, sprintf("%s.mat", name));
disp(matFile);

% check before simulating:
if exist(matFile, 'file')
    fprintf('Found existing results—loading and skipping simulation.\n');
    S = load(matFile, 'simulated_activity_rest_bold_full');
    simulated_activity_rest_bold_full = S.simulated_activity_rest_bold_full;
else
    fprintf('No existing file—running simulation now.\n');

    % % Load 50 fsLR_32k template midthickness surface eigenmodes
    % eigenmodes = dlmread(sprintf('data/examples/fsLR_32k_midthickness-%s_emode_%i.txt', hemisphere, num_modes));
    % eigenvalues = dlmread(sprintf('data/examples/fsLR_32k_midthickness-%s_eval_%i.txt', hemisphere, num_modes));
    
    % Replace above line with the one below and make num_modes = 200 if using the 200 modes provided at data/template_eigenmodes
    eigenmodes = dlmread(sprintf('data/template_eigenmodes/fsLR_32k_midthickness-%s_emode_%i.txt', hemisphere, num_modes));
    eigenvalues = dlmread(sprintf('data/template_eigenmodes/fsLR_32k_midthickness-%s_eval_%i.txt', hemisphere, num_modes));

    % Simulation
    [mode_activity_rest, simulated_activity_rest_full] = ...
        model_neural_waves( ...
            eigenmodes, eigenvalues, ext_input, mu, param, method);
    simulated_activity_rest_full = single(simulated_activity_rest_full);

    [mode_activity_rest, simulated_activity_rest_bold_full] = model_BOLD_balloon(eigenmodes, simulated_activity_rest_full, bal_param, 'Fourier');

    save(matFile, 'simulated_activity_rest_bold_full', '-v7.3');
end
disp("Size of full bold:")
disp(size(simulated_activity_rest_bold_full));


simulated_activity_rest_bold_full = simulated_activity_rest_bold_full(:, end-1200+1:end);
param.T = param.T(end-1200+1:end);
bal_param.T = bal_param.T(end-1200+1:end);

disp("Size of cropped bold:")
disp(size(simulated_activity_rest_bold_full));


analysis_pipeline(simulated_activity_rest_bold_full, name, param.T, cortex, surface_midthickness);

