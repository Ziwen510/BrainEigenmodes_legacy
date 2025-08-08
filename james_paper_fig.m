
addpath(genpath('functions_matlab'));

c_map = load("Rapaeh_color_table_FCD.mat", "c3");
c_map = c_map.c3;

%% MISCELLANEOUS VARIABLES

% NOTE: data provided is only for the below parameters, so please don't change them
hemisphere = 'lh'; 
num_modes = 200;
parc_name = 'Glasser360';

data_empirical_folder = 'data/empirical';
data_results_folder = 'data/results';

%% LOAD SURFACES

% mesh_interest = 'inflated';
% [vertices, faces] = read_vtk(sprintf('data/template_surfaces_volumes/fsLR_32k_%s-%s.vtk', mesh_interest, hemisphere));
% surface_inflated.vertices = vertices';
% surface_inflated.faces = faces';

% mesh_interest = 'midthickness';
% [vertices, faces] = read_vtk(sprintf('data/template_surfaces_volumes/fsLR_32k_%s-%s.vtk', mesh_interest, hemisphere));
% surface_midthickness.vertices = vertices';
% surface_midthickness.faces = faces';

%% LOAD EMPIRICAL AND RESULTS DATA

% cortex and medial wall mask
% cortex = dlmread(sprintf('%s/fsLR_32k_cortex-%s_mask.txt', data_template_surfaces_folder, hemisphere));
% cortex_ind = find(cortex);
% num_vertices = length(cortex);


% RESULTS: model
model_fit = load(sprintf('%s/model_results_Glasser360_%s.mat', data_results_folder, hemisphere), ...
                         'FC_emp', 'FCD_emp', 'FC_model_wave', 'FCD_model_wave', ...
                         'FC_model_mass', 'FCD_model_mass', 'KS');


% disp(size(model_fit.FCD_emp));
% values = model_fit.FCD_emp(:);
%     % Define 100 bins between 0 and 1
% binEdges = linspace(0, 1, 101);  % 101 edges -> 100 bins
% % Compute the PDF using normalized histogram
% pdfArray = histcounts(values, binEdges, 'Normalization', 'probability');
% % Save the PDF array and bin edges to the specified .mat file
% save('./data/empirical/group_FCD_pdf_sliding_Minimal.mat', 'binEdges', 'pdfArray');

FC_emp = load(sprintf('%s/group_FC_ICA+GSR', data_empirical_folder), 'group_fc_mat');
FC_emp = FC_emp.group_fc_mat;
fname = "emp_FC_ICA+GSR";
mask = logical(eye(size(FC_emp)));
FC_emp(mask) = 1;

% FC_emp = model_fit.FC_emp;
% fname = "emp_FC_Minimal";
hFig = figure;  
    imagesc(FC_emp);
    axis square
    colorbar
    clim([-1 1]); 
    colormap(c_map);
    title(fname,'Interpreter','none')
    xlabel('Parcel'), ylabel('Parcel')

    %--- ensure output folder exists & save ------------------------
    savePath = sprintf("./outputs/%s.png", fname);
    [outDir,~,~] = fileparts(savePath);
    if ~isempty(outDir) && ~exist(outDir,'dir')
        mkdir(outDir)
    end
    print(hFig, savePath, '-dpng', '-r300');


%% Plot FCD
% % ut : 694431×1 vector of upper‐triangular FCD entries (i<j)
% % -------------------------------------------------------------------------
% % 1) figure out n such that n*(n-1)/2 == numel(ut)
% ut = model_fit.FCD_model_wave(:, 10);
% n = (1 + sqrt(1 + 8*numel(ut))) / 2;
% if abs(n - round(n)) > 1e-6
%     error('Length of ut (%d) is not n*(n-1)/2 for any integer n.', numel(ut));
% end
% n = round(n);
% 
% % 2) reconstruct the symmetric FCD matrix
% FCD = zeros(n);
% mask = triu(true(n),1);       % logical mask for i<j
% FCD(mask)       = ut;         % fill upper triangle
% FCD = FCD + FCD.';            % mirror to lower triangle
% FCD(1:n+1:end) = 1;           % set diagonal = 1
% 
% % 3) plot
% figure('Color','w','Position',[100 100 600 550]);
% imagesc(FCD);
% axis square ij;               % origin in top‐left
% colormap(turbo);              % perceptually uniform
% colorbar;
% caxis([-1 1]);                % cosine similarity range
% xlabel('Time index');
% ylabel('Time index');
% title('Functional Connectivity Dynamics');
