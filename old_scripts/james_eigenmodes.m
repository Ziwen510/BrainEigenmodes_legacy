%% Load relevant repository MATLAB functions

addpath(genpath('functions_matlab'));

%% MISCELLANEOUS VARIABLES

% NOTE: data provided is only for the below parameters, so please don't change them
hemisphere = 'lh'; 
num_modes = 50;

data_empirical_folder = 'data/empirical';
data_results_folder = 'data/results';
data_figures_folder = 'data/figures_Nature';
data_template_surfaces_folder = 'data/template_surfaces_volumes';
data_template_eigenmodes_folder = 'data/template_eigenmodes';

%% LOAD SURFACES
mesh_interest = 'midthickness';
[vertices, faces] = read_vtk(sprintf('data/template_surfaces_volumes/fsLR_32k_%s-%s.vtk', mesh_interest, hemisphere));
surface_midthickness.vertices = vertices';
surface_midthickness.faces = faces';

%% LOAD EMPIRICAL AND RESULTS DATA


% cortex and medial wall mask
cortex = dlmread(sprintf('%s/fsLR_32k_cortex-%s_mask.txt', data_template_surfaces_folder, hemisphere));
cortex_ind = find(cortex);
num_vertices = length(cortex);

% RESULTS: basis sets     
basis_geometric = dlmread(sprintf('data/examples/fsLR_32k_midthickness-%s_emode_%i.txt', hemisphere, num_modes));



%% DEFINE FIGURE-RELATED PROPERTIES

fontsize_axis = 10;
fontsize_label = 12;


%% EXTENDED DATA FIGURE 1
% bioRxiv Supplementary Figure 1

surface_to_plot = surface_midthickness;

medial_wall = find(~cortex);
mode_list = [1:10,25,50];
num_modes_to_plot = length(mode_list);

factor_x_small = 1.05;
factor_x_big = 2.3;
factor_y = 1.1;
init_x = 0.05;
init_y = 0.01;
length_x = (0.99-init_x)/(factor_x_small + factor_x_big*(4-1) + 1);
length_y = (0.95-init_y)/(factor_y*(num_modes_to_plot-1) + 1);

fig = figure('Position', [200 200 1000 900]);
for basis_ind=1:1
    if basis_ind==1
        data_to_plot = basis_geometric;
        basis_name = 'geometric';
    elseif basis_ind==2
        data_to_plot = basis_connectome;
        basis_name = 'connectome';
    elseif basis_ind==3
        data_to_plot = basis_connectome_density_matched;
        basis_name = {'connectome'; '(density matched)'};
    elseif basis_ind==4
        data_to_plot = basis_EDR;
        basis_name = 'EDR';
    end

    for mode_ind=1:num_modes_to_plot
        mode = mode_list(mode_ind);
        
        data_to_plot(medial_wall,mode) = min(data_to_plot(:,mode))*1.1;
        clims = [min(data_to_plot(:,mode)), max(data_to_plot(:,mode))];
        if clims(2)<=0
            clims(2) = 0.01;
        end

        ax1 = axes('Position', [init_x+factor_x_big*length_x*(basis_ind-1) init_y+factor_y*length_y*(-[mode_ind-num_modes_to_plot]) length_x length_y]);
        obj1 = patch(ax1, 'Vertices', surface_to_plot.vertices, 'Faces', surface_to_plot.faces, 'FaceVertexCData', data_to_plot(:,mode), ...
                   'EdgeColor', 'none', 'FaceColor', 'interp');
        if strcmpi(hemisphere, 'lh')
            view([-90 0]);
        elseif strcmpi(hemisphere, 'rh')
            view([90 0]);
        end
        caxis(clims)
        camlight('headlight')
        material dull
        colormap(ax1,[0.5,0.5,0.5; bluewhitered])
        axis off
        axis image
        
        if basis_ind==1
            annotation(fig, 'textbox', [ax1.Position(1)-0.03, ax1.Position(2)+ax1.Position(4)*0.45, 0.01, 0.01], 'string', {'mode'; sprintf('%i',mode)}, 'edgecolor', 'none', ...
                'fontsize', fontsize_axis, 'fontweight', 'normal', 'horizontalalignment', 'center', 'verticalalignment', 'middle')
        end
        
        ax2 = axes('Position', [init_x+factor_x_small*length_x+factor_x_big*length_x*(basis_ind-1) init_y+factor_y*length_y*(-[mode_ind-num_modes_to_plot]) length_x length_y]);
        obj2 = patch(ax2, 'Vertices', surface_to_plot.vertices, 'Faces', surface_to_plot.faces, 'FaceVertexCData', data_to_plot(:,mode), ...
                   'EdgeColor', 'none', 'FaceColor', 'interp');
        if strcmpi(hemisphere, 'lh')
            view([90 0]);
        elseif strcmpi(hemisphere, 'rh')
            view([-90 0]);
        end
        caxis(clims)
        camlight('headlight')
        material dull
        colormap(ax2,[0.5,0.5,0.5; bluewhitered])
        axis off
        axis image
        
        if mode_ind==1
            annotation(fig, 'textbox', [ax1.Position(1), ax1.Position(2)+ax1.Position(4)*1, ax2.Position(1)-ax1.Position(1)+ax2.Position(3), 0.01], 'string', basis_name, 'edgecolor', 'none', ...
                'fontsize', fontsize_label, 'fontweight', 'normal', 'horizontalalignment', 'center', 'verticalalignment', 'bottom')
        end
    end
end