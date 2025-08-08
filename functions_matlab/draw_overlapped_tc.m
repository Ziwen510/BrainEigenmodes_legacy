function draw_overlapped_tc(ts, fname)
% PLOT_TIMECOURSES  Overlap plot of N time courses and save the figure.
%   plot_timecourses(ts, filename) takes ts of size N×T, plots each of the N
%   time series in ts on a single figure, and saves the plot to filename.

% Get dimensions
[N, T] = size(ts);

% Create figure and plot
figure;
hold on;
for i = 1:N
    plot(1:T, ts(i, :), 'LineWidth', 1);
end
hold off;

% Save figure as PNG at 300 dpi
print(gcf, sprintf("./outputs/%s_overlapped_tc.png", fname), '-dpng', '-r300');
end
