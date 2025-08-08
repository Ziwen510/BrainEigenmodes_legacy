clearvars;         
close all force;
addpath(genpath('functions_matlab'));

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

surface_interest = 'fsLR_32k';
hemisphere = 'lh';
mesh_interest = 'midthickness';

[vertices, faces] = read_vtk(sprintf('data/template_surfaces_volumes/%s_%s-%s.vtk', surface_interest, mesh_interest, hemisphere));
surface_midthickness.vertices = vertices';
surface_midthickness.faces = faces';

% Load cortex mask
cortex = dlmread(sprintf('data/template_surfaces_volumes/%s_cortex-%s_mask.txt', surface_interest, hemisphere));
cortex_ind = find(cortex);

tstep = 180; % in ms

%% Tmax HERE
tmax = 75 * 1000;  % in ms

tspan = [0, tmax];
T = 0:tstep:tmax;
T = T(1:end-1);

%% Simulate resting-state activity

hemisphere = 'lh';
input = 'noise';
Tmax = tmax / 1000;
simulator = 'cortical';
% norm = '_dianorm';
% norm = '_globalnorm';
norm = '';
version = '_debug';
% version = '';

switch simulator
    case 'cortical'
        rate = 1000;
        nu_es = 0.0001 * rate;
        % 
        % nu_ee = 0.00021238 * rate;
        % v0 = 0.00339808 * rate;
        % theta = 0.015 * rate;
        % sigma = 0.0033 * rate;
        
        nu_ee = 0.00006 * rate;
        v0 = 0.0006 * rate;
        theta = 0.0126766 * rate;
        sigma = 0.0038 * rate;
        rmax = '430';
        variable = '_bold';
        boldnorm = '_perframe';
        switch input
            case 'noise'
                seed = 1;
                mean_I = 0;
                ASD = 0.01;
                str = sprintf("white_ASD%.4f_seed%d_mean%.1f", ASD, seed, mean_I);
            case 'pulse'
                cite = 13289;
                % cite = 18245;
                amp = 1000;
                str = sprintf("pulse_%d_I%d", cite, amp);
        end
        S = load_tc(sprintf('my_mats/cortical%s_%s_nuee%.4f_nues%.4f_r30_v3480_v0%.4f_theta%.4f_sigma%.4f_T%.0f_rmax%s%s%s%s.mat', version, str, nu_ee, nu_es, v0, theta, sigma, Tmax, rmax, norm, boldnorm, variable));   % loads only the variable 'tc' into the struct S
    case 'wave'
        K = 255;
        switch input
            case 'noise'
                seed = 1;
                mean_I = 16;
                std = 1;
                str = sprintf("%d_white_std%d_seed%d_mean%.1f", K, std, 1, mean_I);
            case 'pulse'
                cite = 13289;
                % cite = 18245;
                amp = 100;
                mean_I = 0;
                str = sprintf("%d_pulse_%d_I%d_mean%.1f", K, cite, amp, mean_I);
        end
        S = load_tc(sprintf('my_mats/wave%s_r30_v3480_T%.1f.mat', str, Tmax));   % loads only the variable 'tc' into the struct S
end 


simulated_activity_rest = S;
disp(size(simulated_activity_rest));
fname = sprintf("my_%s%s_%s_T%.1f%s", simulator, version, str, Tmax, norm);
disp(fname);


% simulated_actitvity_rest is [Ncortex × T]
% [nc, Tn] = size(simulated_activity_rest);
% Ntotal    = numel(cortex);
% % build full‐size activity: zeros on medial wall
% full_activity = zeros(Ntotal, Tn);
% full_activity(cortex_ind, :) = simulated_activity_rest;
% disp(size(full_activity));

% T_vis = tmax / tstep;
% full_activity = full_activity(:, 1:T_vis);

% %% TC samples before Balloon
% is_cropped = false;
% draw_tc_samples(full_activity, cortex_ind, sprintf("%s_raw_sampled_TC_demeaned", fname), T, 6, is_cropped);
% 
% %% Balloon Model
% full_activity = model_BOLD_balloon_vertex(full_activity, bal_param, 'ODE');
% 
% %% TC samples after Balloon
% % draw_tc_cite(cite, simulated_activity_rest, sprintf("%s_balloon_cite_TC", name), param.T);
% is_cropped = true;
% draw_tc_samples(full_activity, cortex_ind, sprintf("%s_ballon_sampled_TC_demeaned", fname), T, 6, is_cropped);
% 
% %% FC after Balloon
% % draw_corr_matrix(simulated_activity_rest(cortex_ind, :), sprintf('james_raw_phi_FC_modes%d_tstep%d_seed%d_%s_T%d', num_modes, param.tstep, seed, method, param.tmax));
% draw_FC(full_activity, sprintf('%s_balloon_parcellated_FC', fname), hemisphere, is_cropped);



test_pipeline(simulated_activity_rest, cortex_ind, fname, tmax);
% draw_tc_samples(simulated_activity_rest(:,300:1000), cortex_ind, fname, T(300:1000), 10, false);
% test_pipeline_FConly(simulated_activity_rest, cortex_ind, fname);

% =========================================================================
%                      Some visualizations of results                      
% =========================================================================
% 
% draw_tc_cite(15000, full_activity, fname, T(1:T_vis));

%% Snapshot of activity snapshot every 10 ms
% t_interest = [0:10:100];
% t_interest_ind = dsearchn(T', t_interest');
% t_interest_ind = t_interest_ind(1:end-1);
% surface_to_plot = surface_midthickness;
% data_to_plot = full_activity(:, t_interest_ind);
% medial_wall = find(cortex==0);
% with_medial = 0;
% 
% fig = draw_surface_bluewhitered_gallery_dull(surface_to_plot, data_to_plot, hemisphere, medial_wall, with_medial);
% fig.Name = 'Activity snapshots';
% colormap(fig, parula)
% saveas(fig, sprintf("outputs/my_%s.png", str));

%% Video of activity every 0.5 ms (increase this to better see the waves)
% surface_to_plot = surface_midthickness;
% data_to_plot = full_activity;
% t = T;
% tstep_interest = 0.5;
% t_interest = [0:tstep_interest:tmax];
% t_interest = t_interest(10:600);
% is_time_ms = 1;
% medial_wall = find(cortex==0);
% with_medial = 1;
% cmap = parula;
% show_colorbar = 1;
% output_filename = sprintf("outputs/my_%s%s_%s", simulator, version, str);
% save_video = 1;
% 
% fig = video_surface_activity(surface_to_plot, data_to_plot, hemisphere, t, ...
%                              t_interest, is_time_ms, medial_wall, with_medial, ...
%                              cmap, show_colorbar, output_filename, save_video);

