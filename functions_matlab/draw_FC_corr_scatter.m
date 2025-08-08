function corrValue = draw_FC_corr_scatter(sim_FC, fname, emp_version)
    switch emp_version
        case "ICA"
            s2 = load('./data/empirical/group_FC_ICA.mat', 'group_fc_mat');
            y = s2.group_fc_mat;
        case "ICA+GSR"
            s2 = load('./data/empirical/group_FC_ICA+GSR.mat', 'group_fc_mat');
            y = s2.group_fc_mat;
        case "ICA+GSR+WM_CSF_MT_CEN"
            s2 = load('./data/empirical/group_FC_ICA+GSR+WM_CSF_MT_CEN.mat', 'group_fc_mat');
            y = s2.group_fc_mat;
        case "Minimal"
            s2 = load('./data/results/model_results_Glasser360_lh.mat', 'FC_emp');
            y = s2.FC_emp;
    end

    % Only take upper‐triangle entries if these are FC matrices
    if ismatrix(sim_FC) && ~isvector(sim_FC)
        sim_FC = sim_FC(calc_triu_ind(sim_FC));
    end
    if ismatrix(y) && ~isvector(y)
        y = y(calc_triu_ind(y));
    end

    % Fisher‐Z transform
    x = atanh(sim_FC);
    y = atanh(y);

    % Compute Pearson correlation coefficient
    R = corrcoef(x, y);
    corrValue = R(1,2);

    % Create scatter plot
    figure;
    scatter(x, y, 'filled');
    hold on;

    % Fit and plot regression line
    p = polyfit(x, y, 1);
    x_sorted = sort(x);
    y_fit    = polyval(p, x_sorted);
    plot(x_sorted, y_fit, 'r-', 'LineWidth', 2);

    % Labels and grid
    xlabel("Simulated FC", 'Interpreter', 'none');
    ylabel("Empirical FC", 'Interpreter', 'none');
    grid on;
    title(sprintf('%s vs. %s', fname, emp_version), 'Interpreter','none');

    % Add the correlation as text in the upper‐left corner (normalized units)
    txt = sprintf('r = %.2f', corrValue);
    text(0.05, 0.95, txt, 'Units','normalized', ...
         'VerticalAlignment','top', 'FontSize',12, 'BackgroundColor','w', 'EdgeColor','k');

    hold off;

    % Save to PNG
    saveas(gcf, sprintf('./outputs/%s_FC_corr(%s).png', fname, emp_version));
end
