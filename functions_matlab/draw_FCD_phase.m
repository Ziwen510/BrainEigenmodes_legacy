function FCD = draw_FCD_phase(ts, outName)
% COMPUTE_AND_PLOT_FCD  Compute the Functional‑Connectivity Dynamics (FCD),
% plot the FCD matrix, and save it as a PNG.
    c_map = load("Rapaeh_color_table_FCD.mat", "c3");
    c_map = c_map.c3;
    
    ts = ts(:, 22:end);

    % ---- handle defaults ----
    fBand   = [0.04 0.07];
    TR = 0.72;

    [N, T] = size(ts);
    P      = N*(N-1)/2;

    % ---- 1. band‑pass filter (2nd order Butterworth) ----
    fnq = 1/(2*TR);
    [b,a] = butter(2, fBand/fnq);
    % tsF = filtfilt(b, a, detrend(ts','constant'))';  % zero‑phase
    tsF = filtfilt(b, a, ts')';

    % ---- 2. analytic phase via Hilbert ----
    Theta = angle(hilbert(tsF.')).';
    % Theta = angle(hilbert(ts.')).';

    % ---- 3. instantaneous phase‑locking Δ_ij(t) ----
    Delta = zeros(P, T);
    idx   = 0;
    for i = 1:N-1
        for j = i+1:N
            idx = idx + 1;
            Delta(idx,:) = cos(Theta(i,:) - Theta(j,:));
        end
    end

    % ---- 4. normalise each Δ-vector ----
    normDelta = sqrt(sum(Delta.^2, 1));  % 1×T

    % ---- 5. build FCD matrix (cosine similarity) ----
    FCD = (Delta.' * Delta) ./ (normDelta.' * normDelta);
    FCD = max(min(FCD,1),-1);  % clip numerical noise

    % ---- 6. plot ----
    h = figure('Color','w', 'Units','pixels', 'Position', [100 100 600 520]);
    imagesc(FCD);
    axis square ij;
    colormap(turbo);
    colorbar;
    caxis([-1 1]);
    colormap(c_map);
    outName = sprintf("%s_FCD(Phase)", outName);
    title(outName, 'Interpreter','none');
    xlabel('Time index');
    ylabel('Time index');

    % ---- 7. save PNG (300 dpi) ----
    exportgraphics(h, sprintf("./outputs/%s.png", outName), 'Resolution', 300);
    fprintf('Saved FCD figure as ''%s''\n', outName);

end
