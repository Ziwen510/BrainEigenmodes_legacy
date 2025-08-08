function correlation_scatter(matFile1, varName1, label1, matFile2, varName2, label2)
    % Load specified variables from each file
    s1 = load(matFile1, varName1);
    s2 = load(matFile2, varName2);

    % Extract data and ensure as column vectors
    x = s1.(varName1);
    y = s2.(varName2);
    
    % Only take upper‐triangle entries if these are FC matrices
    if ismatrix(x) && ~isvector(x)
        x = x(calc_triu_ind(x));
    end
    if ismatrix(y) && ~isvector(y)
        y = y(calc_triu_ind(y));
    end

    % Fisher‐Z transform
    x = atanh(x);
    y = atanh(y);

    % Check that lengths match
    if numel(x) ~= numel(y)
        error('Error: %s has %d elements, but %s has %d elements. They must be the same length.', ...
              varName1, numel(x), varName2, numel(y));
    end

    % Compute Pearson correlation coefficient
    R = corrcoef(x, y);
    corrValue = R(1,2);

    % Create scatter plot
    figure;
    scatter(x, y, 'filled');
    hold on;

    % Fit and plot red regression line
    p = polyfit(x, y, 1);            % linear fit [slope, intercept]
    x_sorted = sort(x);
    y_fit    = polyval(p, x_sorted);
    plot(x_sorted, y_fit, 'r-', 'LineWidth', 2);

    % Labels, title, grid
    xlabel(label1, 'Interpreter', 'none');
    ylabel(label2, 'Interpreter', 'none');
    % build title string
    titleStr = sprintf('%s vs. %s (r = %.2f)', label1, label2, corrValue);
    title(titleStr, 'Interpreter', 'none');
    grid on;
    hold off;

    % --- Save as PNG with the title as filename (sanitized) ---
    % replace any illegal filename characters with underscore
    safeName = matlab.lang.makeValidName(titleStr);
    % write to PNG in current folder
    saveas(gcf, ['./outputs/' safeName '.png']);
    % alternatively, use print for higher resolution:
    % print(gcf, '-dpng', '-r300', [safeName '.png']);
end



% emp_file = './data/results/model_results_Glasser360_lh.mat';
% emp_var = 'FC_emp';

% emp_file = './data/empirical/group_FC_ICA+GSR.mat';
% emp_var = 'group_fc_mat';

% emp_file = './data/empirical/group_FC_ICA+GSR+WM_CSF_MT_CEN.mat';
% emp_var = 'group_fc_mat';

% my vs. emp
% correlation_scatter('./outputs/my_cortical_debug_white_ASD0.0100_seed1_mean0.0_T865.0_FC.mat', 'FCvec', 'Our Simulated FC', emp_file, emp_var, 'Empirical FC (ICA+GSR+WM_CSF_MT_CEN)');

% james vs. emp
% correlation_scatter('./data/results/model_results_Glasser360_lh.mat', 'FC_model_wave', 'James Simulated FC', emp_file, emp_var, 'Empirical FC (ICA+GSR+WM_CSF_MT_CEN)');

% james vs. my
correlation_scatter('./outputs/james_noise_seed3_mean0_std1.00_Fourier_modes200_T864_tstep0.720_FC.mat', 'FCvec', 'James Simulated FC', './outputs/my_cortical_debug_white_ASD0.0100_seed1_mean0.0_r30_T865.0_FC.mat', 'FCvec', 'Our Simulated FC');

