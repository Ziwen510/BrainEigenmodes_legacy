%% Load relevant repository MATLAB functions

addpath(genpath('functions_matlab'));

%% Load surface files for visualization

surface_interest = 'fsLR_32k';
hemisphere = 'lh';
mesh_interest = 'midthickness';

[vertices, faces] = read_vtk(sprintf('data/template_surfaces_volumes/%s_%s-%s.vtk', surface_interest, mesh_interest, hemisphere));
surface_midthickness.vertices = vertices';
surface_midthickness.faces = faces';

% Load cortex mask
cortex = dlmread(sprintf('data/template_surfaces_volumes/%s_cortex-%s_mask.txt', surface_interest, hemisphere));
cortex_ind = find(cortex);

tstep = 0.1; % in ms
tmax = 100;  % in ms
tspan = [0, tmax];
T = 0:tstep:tmax;
num_modes = 50;

%% Simulate resting-state activity

hemisphere = 'lh';
seed = 3;
S = load(sprintf('my_mats/cortical_white_ASD0.001_seed%d_nuee0.00021238_nues0.0001_mean0.0_r30_v3480_v00.00339808_theta0.015_sigma0.0033_T0.1.mat', seed), 'tc');   % loads only the variable 'tc' into the struct S
simulated_activity_rest = S.tc.';
draw_corr_matrix(simulated_activity_rest, sprintf('my_raw_phi_FC_tstep%d_seed%d', tstep, seed));


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
%     disp(clims);
%     % if clims(2) <= 0
%     %     clims(2) = clims(1) + eps;
%     % end
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
%     clim(clims); colormap([0.5,0.5,0.5; bluewhitered]);
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
%     clim(clims); colormap([0.5,0.5,0.5; bluewhitered]);
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
