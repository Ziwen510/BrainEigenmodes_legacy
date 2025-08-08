function draw_tc_cite(cite, simulated_activity_rest, name, T)
    % build the index range, clipping to valid vertex IDs
    idx_range = (cite-3):(cite+2);
    idx_range = idx_range(idx_range >= 1 & idx_range <= size(simulated_activity_rest,1));
    
    % layout as close to square as possible
    n = numel(idx_range);
    n_cols = 2;
    n_rows = n / n_cols;
    
    fig_raw = figure('Name',name,'NumberTitle','off');
    for k = 1:n
        vi = idx_range(k);
        subplot(n_rows, n_cols, k);
        plot(T, simulated_activity_rest(vi, :), 'LineWidth', 1);
        xlabel('Time (ms)');
        ylabel('Activity');
        title(sprintf('Vertex %d', vi));
        grid on;
    end
    % save to file
    saveas(fig_raw, sprintf("outputs/%s.png", name));
end