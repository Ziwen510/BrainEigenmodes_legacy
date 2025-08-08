function draw_empirical_FCD()

% Load .mat
matfile = "./data/empirical/100206_FCD_phase.mat";
S = load(matfile);
if isfield(S, 'subject')
    subj = S.subject;
else
    subj = 'unknown_subject';
end

% Output directory
outDir = fullfile('.', 'outputs');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

% Define colormap and scaling
c_map = load("Rapaeh_color_table_FCD.mat", "c3");
c_map = c_map.c3;
cmin = -1;
cmax = 1;

% Navigate versions
if ~isfield(S, 'versions')
    error('Loaded file does not contain ''versions'' field.');
end
versions = fieldnames(S.versions);
disp(versions);

for vi = 1:numel(versions)
    version = versions{vi};
    vstruct = S.versions.(version);
    if ~isfield(vstruct, 'per_run')
        warning('Version %s missing per_run; skipping.', version);
        continue;
    end
    runs = fieldnames(vstruct.per_run);
    for ri = 1:numel(runs)
        run = runs{ri};
        runcell = vstruct.per_run.(run);
        % runcell is expected to be a 1x1 struct; extract it safely
        if isempty(runcell) || ~isstruct(runcell)
            warning('Unexpected structure for %s/%s, skipping.', version, run);
            continue;
        end
        runstruct = runcell(1);  % in case it's a 1x1 struct array

        if ~isfield(runstruct, 'FCD_full')
            warning('No FCD_full for %s / %s, skipping.', version, run);
            continue;
        end
        FCD = runstruct.FCD_full;

        if size(FCD,1) ~= size(FCD,2)
            warning('FCD_full for %s/%s is not square; skipping.', version, run);
            continue;
        end

        % Plot
        hF = figure('Visible','off');
        imagesc(FCD);
        axis square;
        colorbar;
        caxis([cmin cmax]);
        colormap(c_map);
        fname = sprintf('%s_%s', version, run);
        title(fname, 'Interpreter', 'none');
        xlabel('Window #');
        ylabel('Window #');

        % Save at 300 dpi
        savePath = fullfile(outDir, sprintf('%s_FCD(Phase).png', fname));
        print(hF, savePath, '-dpng', '-r300');
        close(hF);
        fprintf('Saved %s\n', savePath);
    end
% for vi = 1:numel(versions)
%     version = versions{vi};
%     vstruct = S.versions.(version);
%     run_names = fieldnames(vstruct);
%     for ri = 1:numel(run_names)
%         run = run_names{ri};
%         runstruct_wrapper = vstruct.(run);
%         if isempty(runstruct_wrapper) || ~isstruct(runstruct_wrapper)
%             warning('Skipping %s/%s: not a struct.', version, run);
%             continue;
%         end
%         runstruct = runstruct_wrapper(1);  % unwrap 1x1 struct
% 
%         if ~isfield(runstruct, 'FCD_full')
%             warning('No FCD_full in %s/%s; skipping.', version, run);
%             continue;
%         end
%         FCD = runstruct.FCD_full;
% 
%         if ~ismatrix(FCD) || size(FCD,1) ~= size(FCD,2)
%             warning('FCD_full for %s/%s is not square; skipping.', version, run);
%             continue;
%         end
% 
%         % Plot
%         hF = figure('Visible','off');
%         imagesc(FCD);
%         axis square;
%         colorbar;
%         caxis([cmin cmax]);
%         colormap(c_map);
%         fname = sprintf('%s_%s', version, run);
%         title(fname, 'Interpreter', 'none');
%         xlabel('Window #');
%         ylabel('Window #');
% 
%         % Save at 300 dpi
%         savePath = fullfile(outDir, sprintf('%s_FCD(Phase).png', fname));
%         print(hF, savePath, '-dpng', '-r300');
%         close(hF);
%         fprintf('Saved %s\n', savePath);
%     end
end
end

