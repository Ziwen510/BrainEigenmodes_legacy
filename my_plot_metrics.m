function plot_metrics()
% Load and plot metrics saved in my_mats for r_s = 13:22

% Descriptor for file naming (must match sweep script)
version  = '_debug';
seed     = 1;
mean_I   = 0;
ASD      = 0.01;
str      = sprintf('white_ASD%.4f_seed%d_mean%.1f', ASD, seed, mean_I);
Tmax     = 865.0;  % seconds

% Range of r_s values
rs_vals = 13:22;

% Prepare output directory for figures
output_dir = './outputs';
if ~exist(output_dir,'dir')
    mkdir(output_dir);
end

% Preallocate containers
n = numel(rs_vals);
FC_Minimal             = nan(1,n);
FC_ICA                 = nan(1,n);
FC_ICA_GSR             = nan(1,n);
FC_ICA_GSR_WM          = nan(1,n);
KS_sliding_ICA         = nan(1,n);
KS_sliding_ICA_GSR     = nan(1,n);
KS_sliding_ICA_GSR_WM  = nan(1,n);
KS_phase_Minimal       = nan(1,n);
KS_phase_ICA           = nan(1,n);
KS_phase_ICA_GSR       = nan(1,n);
KS_phase_ICA_GSR_WM    = nan(1,n);

% Loop over r_s
for idx = 1:n
    rs = rs_vals(idx);
    % Construct the base filename (must match metrics_<fname>.mat)
    fname = sprintf('my_cortical%s_%s_r%d_T%.1f', version, str, rs, Tmax);
    metrics_file = fullfile('my_mats', sprintf('metrics_%s.mat', fname));

    if ~isfile(metrics_file)
        warning('Missing metrics file: %s. Skipping r_s=%d', metrics_file, rs);
        continue;
    end
    S = load(metrics_file, 'metrics');
    if ~isfield(S, 'metrics')
        warning('File %s does not contain metrics. Skipping.', metrics_file);
        continue;
    end
    m = S.metrics;

    % Extract FC_corr
    if isfield(m, 'FC_corr')
        fc = m.FC_corr;
        if isfield(fc, 'Minimal'),    FC_Minimal(idx)    = fc.Minimal;    end
        if isfield(fc, 'ICA'),        FC_ICA(idx)        = fc.ICA;        end
        if isfield(fc, 'ICA_GSR'),    FC_ICA_GSR(idx)    = fc.ICA_GSR;    end
        if isfield(fc, 'ICA_GSR_WM_CSF_MT_CEN'), FC_ICA_GSR_WM(idx) = fc.ICA_GSR_WM_CSF_MT_CEN; end
    end

    % Extract KS
    if isfield(m, 'KS')
        ks = m.KS;
        if isfield(ks, 'sliding_ICA'),        KS_sliding_ICA(idx)        = ks.sliding_ICA;        end
        if isfield(ks, 'sliding_ICA_GSR'),    KS_sliding_ICA_GSR(idx)    = ks.sliding_ICA_GSR;    end
        if isfield(ks, 'sliding_ICA_GSR_WM_CSF_MT_CEN'), KS_sliding_ICA_GSR_WM(idx) = ks.sliding_ICA_GSR_WM_CSF_MT_CEN; end
        if isfield(ks, 'phase_Minimal'),      KS_phase_Minimal(idx)      = ks.phase_Minimal;      end
        if isfield(ks, 'phase_ICA'),          KS_phase_ICA(idx)          = ks.phase_ICA;          end
        if isfield(ks, 'phase_ICA_GSR'),      KS_phase_ICA_GSR(idx)      = ks.phase_ICA_GSR;      end
        if isfield(ks, 'phase_ICA_GSR_WM_CSF_MT_CEN'), KS_phase_ICA_GSR_WM(idx) = ks.phase_ICA_GSR_WM_CSF_MT_CEN; end
    end
end

% Plot 1: all FC_corr metrics
figure('Visible','off');
plot(rs_vals, FC_Minimal, '-o', 'DisplayName', 'Minimal'); hold on;
plot(rs_vals, FC_ICA,     '-s', 'DisplayName', 'ICA');
plot(rs_vals, FC_ICA_GSR, '-^', 'DisplayName', 'ICA+GSR');
plot(rs_vals, FC_ICA_GSR_WM, '-d', 'DisplayName', 'ICA+GSR+WM');
xlabel('r_s (mm)'); ylabel('FC correlation');
title('FC Correlation across r_s'); legend('location','best'); grid on;
saveas(gcf, fullfile(output_dir, 'my_FC_corr_vs_rs.png'));
close;

% Plot 2: all KS.sliding metrics
figure('Visible','off');
plot(rs_vals, KS_sliding_ICA, '-o', 'DisplayName', 'sliding\_ICA'); hold on;
plot(rs_vals, KS_sliding_ICA_GSR, '-s', 'DisplayName', 'sliding\_ICA+GSR');
plot(rs_vals, KS_sliding_ICA_GSR_WM, '-^', 'DisplayName', 'sliding\_ICA+GSR+WM');
xlabel('r_s (mm)'); ylabel('KS statistic');
title('Sliding FCD KS across r_s'); legend('location','best'); grid on;
saveas(gcf, fullfile(output_dir, 'my_KS_sliding_vs_rs.png'));
close;

% Plot 3: all KS.phase metrics
figure('Visible','off');
plot(rs_vals, KS_phase_Minimal, '-o', 'DisplayName', 'phase\_Minimal'); hold on;
plot(rs_vals, KS_phase_ICA,     '-s', 'DisplayName', 'phase\_ICA');
plot(rs_vals, KS_phase_ICA_GSR, '-^', 'DisplayName', 'phase\_ICA+GSR');
plot(rs_vals, KS_phase_ICA_GSR_WM, '-d', 'DisplayName', 'phase\_ICA+GSR+WM');
xlabel('r_s (mm)'); ylabel('KS statistic');
title('Phase FCD KS across r_s'); legend('location','best'); grid on;
saveas(gcf, fullfile(output_dir, 'my_KS_phase_vs_rs.png'));
close;

% Plot 4: FC_corr.ICA & KS.sliding_ICA
figure('Visible','off');
yyaxis left;  plot(rs_vals, FC_ICA, '-o', 'DisplayName', 'FC (ICA)'); ylabel('FC (ICA)');
yyaxis right; plot(rs_vals, KS_sliding_ICA, '-s', 'DisplayName', 'KS sliding (ICA)'); ylabel('KS sliding (ICA)');
xlabel('r_s (mm)'); title('FC Correlation & Sliding KS (ICA)'); legend('location','best'); grid on;
saveas(gcf, fullfile(output_dir, 'my_FC_vs_KS_sliding_ICA.png'));
close;

% Plot 5: FC_corr.Minimal & KS.phase_Minimal
figure('Visible','off');
yyaxis left;  plot(rs_vals, FC_Minimal, '-o', 'DisplayName', 'FC (Minimal)'); ylabel('FC (Minimal)');
yyaxis right; plot(rs_vals, KS_phase_Minimal, '-s', 'DisplayName', 'KS phase (Minimal)'); ylabel('KS phase (Minimal)');
xlabel('r_s (mm)'); title('FC & Phase KS (Minimal)'); legend('location','best'); grid on;
saveas(gcf, fullfile(output_dir, 'my_FC_Minimal_vs_KS_phase_Minimal.png'));
close;

% Plot 6: KS_sliding_ICA, KS_phase_ICA & KS_phase_Minimal
figure('Visible','off');
plot(rs_vals, KS_sliding_ICA, '-o', 'DisplayName', 'KS sliding (ICA)'); hold on;
plot(rs_vals, KS_phase_ICA,   '-s', 'DisplayName', 'KS phase (ICA)');
plot(rs_vals, KS_phase_Minimal,'-d', 'DisplayName', 'KS phase (Minimal)');
xlabel('r_s (mm)'); ylabel('KS statistic');
title('KS summary vs r_s'); legend('location','best'); grid on;
saveas(gcf, fullfile(output_dir, 'my_KS_summary_vs_rs.png'));
close;
end
