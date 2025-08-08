function plot_metrics()
% Assumes files named like: metrics_james_noise_rs<r_s>_Fourier_modes200_T864_tstep0.720.mat
% Adjust `input`, `method`, `num_modes`, `tmax`, `tstep` if you used different naming.

input = 'noise';
method = 'Fourier';
num_modes = 200;
tmax = 864;
tstep = 0.72;
source = "my";

rs_vals = 10:60;

output_dir = './outputs';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Preallocate containers
n = numel(rs_vals);
FC_Minimal = nan(1,n);
FC_ICA = nan(1,n);
FC_ICA_GSR = nan(1,n);
FC_ICA_GSR_WM = nan(1,n);

KS_sliding_ICA = nan(1,n);
KS_sliding_ICA_GSR = nan(1,n);
KS_sliding_ICA_GSR_WM = nan(1,n);

KS_phase_Minimal = nan(1,n);
KS_phase_ICA = nan(1,n);
KS_phase_ICA_GSR = nan(1,n);
KS_phase_ICA_GSR_WM = nan(1,n);

for idx = 1:n
    rs = rs_vals(idx);
    name_rs = sprintf("james_%s_rs%d_%s_modes%d_T%d_tstep%.3f", input, rs, method, num_modes, tmax, tstep);
    metrics_fname = sprintf("./james_mats/metrics_%s.mat", name_rs);
    if ~isfile(metrics_fname)
        warning("Missing metrics file: %s. Skipping r_s=%d", metrics_fname, rs);
        continue;
    end
    S = load(metrics_fname, 'metrics');
    if ~isfield(S, 'metrics')
        warning("File %s does not contain 'metrics' variable. Skipping.", metrics_fname);
        continue;
    end
    m = S.metrics;

    % FC_corr
    if isfield(m, 'FC_corr')
        fc = m.FC_corr;
        if isfield(fc, 'Minimal'), FC_Minimal(idx) = fc.Minimal; end
        if isfield(fc, 'ICA'), FC_ICA(idx) = fc.ICA; end
        if isfield(fc, 'ICA_GSR'), FC_ICA_GSR(idx) = fc.ICA_GSR; end
        if isfield(fc, 'ICA_GSR_WM_CSF_MT_CEN'), FC_ICA_GSR_WM(idx) = fc.ICA_GSR_WM_CSF_MT_CEN; end
    end

    % KS
    if isfield(m, 'KS')
        ks = m.KS;
        if isfield(ks, 'sliding_ICA'), KS_sliding_ICA(idx) = ks.sliding_ICA; end
        if isfield(ks, 'sliding_ICA_GSR'), KS_sliding_ICA_GSR(idx) = ks.sliding_ICA_GSR; end
        if isfield(ks, 'sliding_ICA_GSR_WM_CSF_MT_CEN'), KS_sliding_ICA_GSR_WM(idx) = ks.sliding_ICA_GSR_WM_CSF_MT_CEN; end

        if isfield(ks, 'phase_Minimal'), KS_phase_Minimal(idx) = ks.phase_Minimal; end
        if isfield(ks, 'phase_ICA'), KS_phase_ICA(idx) = ks.phase_ICA; end
        if isfield(ks, 'phase_ICA_GSR'), KS_phase_ICA_GSR(idx) = ks.phase_ICA_GSR; end
        if isfield(ks, 'phase_ICA_GSR_WM_CSF_MT_CEN'), KS_phase_ICA_GSR_WM(idx) = ks.phase_ICA_GSR_WM_CSF_MT_CEN; end
    end
end

% Plot 1: all FC_corr metrics
figure;
plot(rs_vals, FC_Minimal, '-o', 'DisplayName', 'Minimal'); hold on;
plot(rs_vals, FC_ICA, '-s', 'DisplayName', 'ICA');
plot(rs_vals, FC_ICA_GSR, '-^', 'DisplayName', 'ICA+GSR');
plot(rs_vals, FC_ICA_GSR_WM, '-d', 'DisplayName', 'ICA+GSR+Others');
xlabel('r_s (mm)');
ylabel('FC correlation');
title('FC Correlation across r_s');
legend('location','best');
grid on;
saveas(gcf, fullfile(output_dir, 'FC_corr_vs_rs.png'));
close;

% Plot 2: all KS.sliding metrics
figure;
plot(rs_vals, KS_sliding_ICA, '-o', 'DisplayName', 'sliding\_ICA'); hold on;
% plot(rs_vals, KS_sliding_ICA_GSR, '-s', 'DisplayName', 'sliding\_ICA+GSR');
% plot(rs_vals, KS_sliding_ICA_GSR_WM, '-^', 'DisplayName', 'sliding\_ICA+GSR+Others');
xlabel('r_s (mm)');
ylabel('KS statistic');
title('Sliding FCD KS across r_s');
legend('location','best');
grid on;
saveas(gcf, fullfile(output_dir, 'KS_sliding_vs_rs.png'));
close;

% Plot 3: all KS.phase metrics
figure;
plot(rs_vals, KS_phase_Minimal, '-o', 'DisplayName', 'phase\_Minimal'); hold on;
plot(rs_vals, KS_phase_ICA, '-s', 'DisplayName', 'phase\_ICA');
% plot(rs_vals, KS_phase_ICA_GSR, '-^', 'DisplayName', 'phase\_ICA+GSR');
% plot(rs_vals, KS_phase_ICA_GSR_WM, '-d', 'DisplayName', 'phase\_ICA+GSR+Others');
xlabel('r_s (mm)');
ylabel('KS statistic');
title('Phase FCD KS across r_s');
legend('location','best');
grid on;
saveas(gcf, fullfile(output_dir, 'KS_phase_vs_rs.png'));
close;

% Plot 4: FC_corr.ICA & KS.sliding_ICA
figure;
yyaxis left
plot(rs_vals, FC_ICA, '-o', 'DisplayName', 'FC Correlation (ICA)');
ylabel('FC Correlation (ICA)');
yyaxis right
plot(rs_vals, KS_sliding_ICA, '-s', 'DisplayName', 'Sliding KS (ICA)');
ylabel('Sliding KS (ICA)');
xlabel('r_s (mm)');
title('FC Correlation & Sliding KS (ICA) across r_s');
legend('location','best');
grid on;
saveas(gcf, fullfile(output_dir, 'FC_vs_KS_sliding_ICA.png'));
close;

% Plot 5: FC_corr.Minimal & KS.phase_Minimal
figure;
yyaxis left
plot(rs_vals, FC_Minimal, '-o', 'DisplayName', 'FC Correlation (Minimal)');
ylabel('FC Correlation (Minimal)');
yyaxis right
plot(rs_vals, KS_phase_Minimal, '-s', 'DisplayName', 'Phase KS (Minimal)');
ylabel('Phase KS (Minimal)');
xlabel('r_s (mm)');
title('FC Correlation & Phase KS (Minimal) across r_s');
legend('location','best');
grid on;
saveas(gcf, fullfile(output_dir, 'FC_Minimal_vs_KS_phase_Minimal.png'));
close;

% Plot 6: KS_sliding_ICA, KS_phase_ICA & KS_phase_Minimal
figure('Visible','off');
plot(rs_vals, KS_sliding_ICA, '-o', 'DisplayName', 'Sliding KS (ICA)'); hold on;
plot(rs_vals, KS_phase_ICA, '-s', 'DisplayName', 'Phase KS (ICA)');
plot(rs_vals, KS_phase_Minimal, '-d', 'DisplayName', 'Phase KS (Minimal)');
xlabel('r_s (mm)');
ylabel('KS statistic');
title('KS with sliding\_ICA, phase\_ICA & phase\_Minimal vs r_s');
legend('location','best');
grid on;
saveas(gcf, fullfile(output_dir, 'KS_display.png'));
close;

end
