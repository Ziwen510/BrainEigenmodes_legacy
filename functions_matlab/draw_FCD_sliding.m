function FCD = draw_FCD_sliding(data_parc, window_size, fname)
% input: N*T
c_map = load("Rapaeh_color_table_FCD.mat", "c3");
c_map = c_map.c3;

%% compute number of windows
[~, T] = size(data_parc);
tlen   = T - window_size + 1;
disp("FCD tlen: ");
disp(tlen);
if tlen < 2
    error('window_size must be < number of timepoints.');
end

%% preallocate feature matrix (windows × FC edges)
nROI      = size(data_parc,1);
fc_edge   = nROI*(nROI-1)/2;
fc_vecs   = zeros(tlen, fc_edge);
disp("nROI: ");
disp(nROI);

%% upper‐triangle mask for vectorization
mask_ut = triu(true(nROI),1);

%% sliding‐window FC
for w = 1:tlen
    seg = data_parc(:, w:(w+window_size-1));   % [nROI×window_size]
    C   = corrcoef(seg');                       % [nROI×nROI]
    fc_vecs(w,:) = C(mask_ut)';                 % 1×fc_edge
end

%% compute FCD matrix (correlation of FC vectors)
FCD = corrcoef(fc_vecs');  % [tlen×tlen]
disp("shape of Sliding FCD matrix");
disp(size(FCD));

%% plot & save
hF = figure();
imagesc(FCD);
axis square;
colorbar;
caxis([0.1 1]);
colormap(c_map);
title(fname,'Interpreter','none');
xlabel('Window #');
ylabel('Window #');

% ensure output folder
outDir = fullfile('.', 'outputs');
if ~exist(outDir,'dir')
    mkdir(outDir);
end

% save at 300 dpi
savePath = fullfile(outDir, sprintf('%s_FCD(Sliding).png', fname));
print(hF, savePath, '-dpng', '-r300');


end
