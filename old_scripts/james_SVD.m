% demo_wave_model_svd_visualization.m
% ------------------------------------------------------------
% This script simulates resting-state wave activity on a cortical
% surface, performs SVD on the simulated time course, and then
% visualizes the first 50 spatial modes (SVD basis) exactly as in
% the Extended-Data-Fig.1 “geometric” gallery.
%
% Original demo: James Pang, Monash University, 2022
% Adapted for SVD visualization: Wang Ziwen, 2025
% ------------------------------------------------------------

%% 1. Add paths and define parameters
addpath(genpath('functions_matlab'));

surface_interest          = 'fsLR_32k';
hemisphere                = 'lh';
mesh_interest             = 'midthickness';
num_modes                 = 50;         % number of modes to simulate & recover

%% 2. Load surface geometry
[verts, faces] = read_vtk( ...
    sprintf('data/template_surfaces_volumes/%s_%s-%s.vtk', ...
    surface_interest, mesh_interest, hemisphere) );
surface_midthickness.vertices = verts';
surface_midthickness.faces    = faces';

%% 3. Load cortex mask
cortex    = dlmread( ...
    sprintf('data/template_surfaces_volumes/%s_cortex-%s_mask.txt', ...
    surface_interest, hemisphere) );
cortex_ind = find(cortex);

%% 4. Simulate resting-state activity (50 modes)
eigenmodes  = dlmread( ...
    sprintf('data/examples/fsLR_32k_midthickness-%s_emode_%i.txt', ...
    hemisphere, num_modes) );
eigenvalues = dlmread( ...
    sprintf('data/examples/fsLR_32k_midthickness-%s_eval_%i.txt', ...
    hemisphere, num_modes) );

param         = loadParameters_wave_func;
param.tstep   = 0.1;            % ms
param.tmax    = 5000;            % ms
param.T       = 0:param.tstep:param.tmax;
param.is_time_ms = 1;
param.r_s     = 30;             % mm
param.gamma_s = 116 * 1e-3;     % s^-1 → ms^-1

method = 'Fourier';  % use 'Fourier' for long T

seed = 1; 
rng(seed);  % reproducible white noise
ext_input = randn(size(eigenmodes,1), numel(param.T));

[mode_activity_rest, simulated_activity_rest] = ...
    model_neural_waves(eigenmodes, eigenvalues, ext_input, param, method);

draw_corr_matrix(simulated_activity_rest(cortex_ind, :), sprintf('james_raw_phi_FC_modes%d_tstep%d_seed%d_%s_T%d', num_modes, param.tstep, seed, method, param.tmax));

% %% 5. Perform SVD on the full V×T data
% % simulated_activity_rest is V×T
% [U_svd, S_svd, V_svd] = svd(simulated_activity_rest, 'econ');
% basis_svd = U_svd(:, 1:num_modes);  % recover first 50 modes
% 
% %% 6. Visualize the recovered SVD modes
% surface_to_plot    = surface_midthickness;
% medial_wall        = find(~cortex);      % indices to mask out
% mode_list          = [1:10, 25, 50];
% num_modes_to_plot  = numel(mode_list);
% 
% % Layout parameters (same logic as Extended-Data-Fig1)
% factor_x_small = 1.05;
% factor_x_big   = 2.3;
% factor_y       = 1.1;
% init_x         = 0.05;
% init_y         = 0.01;
% length_x = (0.99 - init_x) / (factor_x_small + 1);
% length_y = (0.95 - init_y) / (factor_y*(num_modes_to_plot-1) + 1);
% 
% fig = figure('Position',[200 200 1000 900]);
% for mi = 1:num_modes_to_plot
%     mode = mode_list(mi);
% 
%     % extract and mask the spatial map
%     d = basis_svd(:, mode);
%     d(medial_wall) = min(d)*1.1;
% 
%     clims = [min(d), max(d)];
%     if clims(2) <= 0
%         clims(2) = clims(1) + eps;
%     end
% 
%     % — Left (lateral) view
%     ax1 = axes('Position',[
%         init_x, ...
%         init_y + (num_modes_to_plot-mi)*factor_y*length_y, ...
%         length_x, length_y]);
%     patch('Vertices', surface_to_plot.vertices, ...
%           'Faces',    surface_to_plot.faces, ...
%           'FaceVertexCData', d, ...
%           'EdgeColor','none','FaceColor','interp');
%     view([-90 0]); camlight headlight; material dull;
%     caxis(clims); colormap([0.5,0.5,0.5; bluewhitered]);
%     axis off; axis image;
% 
%     % — Right (medial) view
%     ax2 = axes('Position',[
%         init_x + factor_x_small*length_x, ...
%         init_y + (num_modes_to_plot-mi)*factor_y*length_y, ...
%         length_x, length_y]);
%     patch('Vertices', surface_to_plot.vertices, ...
%           'Faces',    surface_to_plot.faces, ...
%           'FaceVertexCData', d, ...
%           'EdgeColor','none','FaceColor','interp');
%     view([90 0]); camlight headlight; material dull;
%     caxis(clims); colormap([0.5,0.5,0.5; bluewhitered]);
%     axis off; axis image;
% 
%     % Mode number annotation
%     annotation(fig,'textbox',[
%         ax1.Position(1)-0.03, ...
%         ax1.Position(2)+ax1.Position(4)/2, ...
%         0.01, 0.01], ...
%         'String',sprintf('mode %d',mode), ...
%         'EdgeColor','none','FontSize',10, ...
%         'HorizontalAlignment','center');
% 
%     % Basis name at top of first row
%     if mi==1
%         annotation(fig,'textbox',[
%             ax1.Position(1), ...
%             ax1.Position(2)+ax1.Position(4)*1.02, ...
%             (ax2.Position(1)+ax2.Position(3)-ax1.Position(1)), ...
%             0.02], ...
%             'String','SVD modes', ...
%             'EdgeColor','none','FontSize',12, ...
%             'HorizontalAlignment','center', ...
%             'VerticalAlignment','bottom');
%     end
% end
