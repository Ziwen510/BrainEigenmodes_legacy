clearvars;
close all force;
addpath(genpath('functions_matlab'));
set(0, 'DefaultFigureVisible', 'off');

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
%         Simulation sweep over r_s   
% =========================================================================

tmax = 864; % in s
tstep = 0.72;

input = 'noise';  % or 'pulse'
method = 'Fourier';

% build common output directory
outputDir = fullfile('.', 'james_mats');
if ~exist(outputDir, 'dir'), mkdir(outputDir); end

% initialize result container
rs_vals = 10:60;
results = struct('r_s', {}, 'FC_corr', {}, 'KS', {});

% fixed external input (does not depend on r_s)
% prepare ext_input once
switch input
    case 'noise'
        seed = 1;
        rng(seed);
        mu = 0;
        stdv = 1;
        ext_input = randn(length(cortex), ceil(tmax / tstep) + 1) * stdv + mu; % full length; will be cropped later
        str = sprintf("seed%d_mean%d_std%.2f", seed, mu, stdv);
    case 'pulse'
        cite = 15000;
        amp = 100;
        mu = 0;
        ext_input = ones(length(cortex), ceil(tmax / tstep) + 1) * mu;
        ext_input(cite, 1:10) = ext_input(cite, 1:10) + amp;
        str = sprintf("cite%d_amp%d_mean%d", cite, amp, mu);
end

for idx = 1:numel(rs_vals)
    rs = rs_vals(idx);
    fprintf("===== Running / loading for r_s = %d =====\n", rs);

    % Reload parameters fresh each iteration to avoid carryover
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

    param.r_s = rs;      % variable of interest
    param.gamma_s = 116; % (default) in s^-1
    if param.is_time_ms == 1
        param.gamma_s = 116 * 1e-3;
    end

    name_rs = sprintf("james_%s_rs%d_%s_modes%d_T%d_tstep%.3f", input, rs, method, num_modes, param.tmax, param.tstep);
    matFile_rs = fullfile(outputDir, sprintf("%s.mat", name_rs));
    disp(matFile_rs);

    if exist(matFile_rs, 'file')
        fprintf('Found existing results—loading and skipping simulation for r_s=%d.\n', rs);
        S = load(matFile_rs, 'simulated_activity_rest_bold_full');
        simulated_activity_rest_bold_full = S.simulated_activity_rest_bold_full;
    else
        fprintf('No existing file—running simulation now for r_s=%d.\n', rs);
        % load eigenmodes / eigenvalues
        eigenmodes = dlmread(sprintf('data/template_eigenmodes/fsLR_32k_midthickness-%s_emode_%i.txt', hemisphere, num_modes));
        eigenvalues = dlmread(sprintf('data/template_eigenmodes/fsLR_32k_midthickness-%s_eval_%i.txt', hemisphere, num_modes));

        % Simulation
        [mode_activity_rest, simulated_activity_rest_full] = ...
            model_neural_waves( ...
                eigenmodes, eigenvalues, ext_input, mu, param, method);
        simulated_activity_rest_full = single(simulated_activity_rest_full);

        [mode_activity_rest, simulated_activity_rest_bold_full] = model_BOLD_balloon(eigenmodes, simulated_activity_rest_full, bal_param, 'Fourier');

        save(matFile_rs, 'simulated_activity_rest_bold_full', '-v7.3');
    end

    % Crop last 1200 time points
    simulated_activity_rest_bold_full = simulated_activity_rest_bold_full(:, end-1200+1:end);
    param.T = param.T(end-1200+1:end);
    bal_param.T = bal_param.T(end-1200+1:end);

    % Run modified pipeline that returns metrics
    metrics = analysis_pipeline(simulated_activity_rest_bold_full, name_rs, param.T, cortex, surface_midthickness);

    % attach r_s for bookkeeping
    metrics.r_s = rs;
    results(end+1) = metrics;

    % optionally save per-rs metrics as well
    sweep_matname = fullfile(outputDir, sprintf("metrics_%s.mat", name_rs));
    save(sweep_matname, 'metrics', '-v7.3');
end

% Save aggregated sweep
aggregate_name = fullfile(outputDir, sprintf("james_%s_%s_modes%d_rsweep.mat", input, method, num_modes));
save(aggregate_name, 'results', 'rs_vals', '-v7.3');
fprintf("Saved aggregated results to %s\n", aggregate_name);
