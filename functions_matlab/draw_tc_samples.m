function draw_tc_samples(simulated_activity_rest, cortex_ind, name, T, num_samples, is_cropped)
    % draw_tc_samples  Plot time‐courses of evenly spaced cortex vertices
    %
    %   draw_tc_samples(simulated_activity_rest, cortex_ind, name, T, num_samples)
    %
    %   Inputs:
    %     simulated_activity_rest  [N×M] matrix of N vertices × M time‐points
    %     cortex_ind               vector of valid vertex indices (subset of 1…N)
    %     name                     string used to name the output file
    %     T                        vector of time‐points (length M)
    %     num_samples              number of vertices to sample & plot
    %
    %   Example:
    %     draw_tc_samples(activity, cortex_ind, 'run1', time_vector, 8);
    % disp(T);
    if nargin < 5 || isempty(num_samples)
        num_samples = 6;   % default to 6 samples
    end

    % ensure cortex_ind is sorted and within valid range
    cortex_ind = cortex_ind(:);
    cortex_ind = cortex_ind(cortex_ind >= 1 & cortex_ind <= size(simulated_activity_rest,1));

    % how many cortex vertices we have
    Nc = numel(cortex_ind);
    if num_samples > Nc
        warning('Requested %d samples but only %d cortex vertices available. Reducing to %d.', ...
                num_samples, Nc, Nc);
        num_samples = Nc;
    end

    % pick indices evenly spaced within cortex_ind
    sample_positions = round(linspace(1, 20000, num_samples));
    idx_samples = cortex_ind(sample_positions);

    % determine subplot grid (near square)
    n = numel(idx_samples);
    n_cols = 2;
    n_rows = ceil(n / n_cols);

    % create figure
    fig_raw = figure('Name', name, ...
                     'NumberTitle', 'off');

    if is_cropped
        T = T(:, 300:end);
        simulated_activity_rest = simulated_activity_rest(:, 300:end);
        name = sprintf("%s_cropped", name);
    end
    % plot each sampled cortex vertex
    for k = 1:n
        vi = idx_samples(k);
        subplot(n_rows, n_cols, k);
        plot(T, simulated_activity_rest(vi, :), 'LineWidth', 1);
        xlim([T(1), T(end)]);
        xlabel('Time (s)');
        ylabel('Activity');
        title(sprintf('Vertex %d', vi));
        grid on;
    end

    % save to file (as PNG)
    saveas(fig_raw, fullfile('outputs', sprintf('%s.png', name)));
end
