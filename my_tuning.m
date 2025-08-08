%% Sweep r_s from 13 to 22 and save metrics per r_s
clearvars;
close all force;
addpath(genpath('functions_matlab'));
set(0,'DefaultFigureVisible','off');


function tc = load_tc(filename)
%LOAD_TC   Load the time-course matrix 'tc' from a MAT-file.
%   tc = LOAD_TC(filename)  
%   - First attempts S = load(filename,'tc'); tc = S.tc;
%   - If that errors (e.g. data >4GB in v7 MAT), it will try h5read(filename,'/tc').
%
%   Example:
%     tc = load_tc('my_mats/cortical_str_nuee0.2000_nues0.1000.mat');

    % validate input
    if ~isfile(filename)
        error('File not found: %s', filename);
    end

    % try standard MAT-file load
    try
        S = load(filename, 'tc');
        tc = S.tc';
        return
    catch ME
        warning('Standard load failed: %s\nFalling back to HDF5…', ME.message);
    end

    % fallback to HDF5 read
    try
        tc = h5read(filename, '/tc');
    catch H5E
        error('HDF5 read failed for %s:\n%s', filename, H5E.message);
    end
end

% Build output directory
outputDir = fullfile('.', 'my_mats');
if ~exist(outputDir,'dir')
    mkdir(outputDir);
end

% Load surface geometry and cortex mask
surface_interest = 'fsLR_32k';
hemisphere = 'lh';
mesh_interest = 'midthickness';

cortex = dlmread(sprintf('data/template_surfaces_volumes/%s_cortex-%s_mask.txt', surface_interest, hemisphere));
cortex_ind = find(cortex);
[vertices, faces] = read_vtk(sprintf('data/template_surfaces_volumes/%s_%s-%s.vtk', surface_interest, mesh_interest, hemisphere));
surface_midthickness.vertices = vertices';
surface_midthickness.faces    = faces';

% Simulation timing parameters
tmax    = 865 * 1000;   % ms
tstep   = 720;          % ms
Tmax    = tmax / 1000;  % s for naming
warmup  = 1440;         % ms
crop_dt = 720;          % ms

% Model parameters (constant across r_s)
input    = 'noise';
rate     = 1000;
nu_es    = 0.0001 * rate;
nu_ee    = 0.00006 * rate;
v0       = 0.0006 * rate;
theta    = 0.0126766 * rate;
sigma    = 0.0038 * rate;
version  = '_debug';
norm     = '';
boldnorm = '';
variable = '_bold';
rmax     = '430';
cropping = '';

% Prepare input descriptor string
switch input
    case 'noise'
        seed   = 1;
        mean_I = 0;
        ASD    = 0.01;
        str    = sprintf('white_ASD%.4f_seed%d_mean%.1f', ASD, seed, mean_I);
    case 'pulse'
        cite = 13289;
        amp  = 1000;
        str  = sprintf('pulse_%d_I%d', cite, amp);
end

% Time vector for diagnostic
T = 0:tstep:tmax;
disp('Original length of T:'), disp(numel(T));

% r_s sweep values
rs_vals = 13:22;
results = struct('r_s', {}, 'FC_corr', {}, 'KS', {});

for i = 1:numel(rs_vals)
    r_s = rs_vals(i);
    v   = r_s * 116;
    
    % Construct base filename/fname
    fname = sprintf('my_cortical%s_%s_r%d_T%.1f%s%s', version, str, r_s, Tmax, norm, cropping);
    fprintf('--- Processing r_s = %d (%s) ---\n', r_s, fname);
    
    % Load time course using helper
    tc_file = fullfile(outputDir, sprintf(...
        'cortical%s_%s_nuee%.4f_nues%.4f_r%d_v%d_v0%.4f_theta%.4f_sigma%.4f_T%.0f_rmax%s%s%s%s.mat', ...
        version, str, nu_ee, nu_es, r_s, v, v0, theta, sigma, Tmax, rmax, norm, boldnorm, variable));
    tc = load_tc(tc_file);
    
    % Crop and downsample
    switch cropping
        case '_cropped'
            T_target = warmup:crop_dt:217*1000;
        otherwise
            T_target = warmup:crop_dt:tmax;
    end
    idx_T   = round(T_target./tstep) + 1;
    tc_down = tc(:, idx_T);   % [Ncortex x M]
    
    % Expand to full cortex grid
    [~, M]     = size(tc_down);
    full_tc    = zeros(numel(cortex), M);
    full_tc(cortex_ind, :) = tc_down;
    
    % Run analysis pipeline to get metrics struct
    metrics         = analysis_pipeline(full_tc, fname, T_target, cortex, surface_midthickness);
    metrics.r_s     = r_s;
    results(end+1) = metrics;
    
    % Save per-r_s metrics
    out_metrics_fn = fullfile(outputDir, sprintf('metrics_%s.mat', fname));
    save(out_metrics_fn, 'metrics', '-v7.3');
end

% Save aggregated results
agg_fn = fullfile(outputDir, 'metrics_rsweep.mat');
save(agg_fn, 'results', 'rs_vals', '-v7.3');
fprintf('Saved aggregated metrics to %s\n', agg_fn);
