function C = draw_corr_matrix(timecourse, name)
    % validate input
    if ndims(timecourse)~=2
        error('Input must be a 2D matrix of size N×T.');
    end
    disp(size(timecourse));

    % Step 1: demean per frame (column)
    % compute mean across rows for each time point
    frameMeans = mean(timecourse, 1);        % 1×T
    % subtract from each row
    X = timecourse - frameMeans;             % N×T

    % Step 2: compute correlation matrix across rows
    % corrcoef operates on columns, so transpose
    C = corrcoef(X');                         % N×N

    % Step 3: plot
    h = figure('Visible','off');              % create invisible figure
    imagesc(C);
    axis square tight;
    colorbar;
    title(name,'Interpreter','none');
    xlabel('Vertex index');
    ylabel('Vertex index');
    colormap(parula);                         % or choose your favorite

    % Step 4: save
    fname = sprintf('./outputs/%s.png', name);
    % you can also use saveas(h, fname) or print:
    saveas(h, fname);

    % close if you don't want it open in MATLAB
    close(h);

    fprintf('Correlation matrix computed and plot saved to %s\n', fname);
end
