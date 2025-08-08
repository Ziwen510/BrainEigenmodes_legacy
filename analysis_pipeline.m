function metrics = analysis_pipeline(bold_full, name, T, cortex, surface)
    cortex_ind = find(cortex);
    medial_wall = find(~cortex);

    metrics = struct();

    %% Mean BOLD Plot
    mean_activity = mean(bold_full, 1);
    h = figure;
    plot(mean_activity, 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel('Mean BOLD');
    title('Mean BOLD Across Vertices');
    grid on;
    print(h, sprintf('./outputs/%s_mean_BOLD.png', name), '-dpng', '-r300');
    close(h);

    %% Sample Time Courses
    draw_tc_samples(bold_full, cortex_ind, sprintf("%s_BOLD_TC_samples", name), T, 8, false);

    %% Parcellate
    data_parc = my_parcellate(bold_full, 'lh');
    disp("Size of parcelled bold:")
    disp(size(data_parc));

    % draw_overlapped_tc(data_parc, name);

    %% FC and correlations
    FC_vec = draw_FC(data_parc, name);

    % initialize container for correlations
    FC_corr = struct();
    FC_corr.Minimal = draw_FC_corr_scatter(FC_vec, name, 'Minimal');
    FC_corr.ICA = draw_FC_corr_scatter(FC_vec, name, 'ICA');
    FC_corr.ICA_GSR = draw_FC_corr_scatter(FC_vec, name, 'ICA+GSR');
    FC_corr.ICA_GSR_WM_CSF_MT_CEN = draw_FC_corr_scatter(FC_vec, name, 'ICA+GSR+WM_CSF_MT_CEN');

    %% FCD + KS statistics
    KS = struct();

    % sliding, ICA variants
    FCD_sliding = draw_FCD_sliding(data_parc, 83, name);
    KS.sliding_ICA = draw_FCD_pdf_KS(FCD_sliding, "sliding_ICA", name);
    % KS.sliding_ICA_GSR = draw_FCD_pdf_KS(FCD_sliding, "sliding_ICA+GSR", name);
    % KS.sliding_ICA_GSR_WM_CSF_MT_CEN = draw_FCD_pdf_KS(FCD_sliding, "sliding_ICA+GSR+WM_CSF_MT_CEN", name);

    % phase variants
    FCD_phase = draw_FCD_phase(data_parc, name);
    KS.phase_Minimal = draw_FCD_pdf_KS(FCD_phase, "phase_Minimal", name);
    KS.phase_ICA = draw_FCD_pdf_KS(FCD_phase, "phase_ICA", name);
    % KS.phase_ICA_GSR = draw_FCD_pdf_KS(FCD_phase, "phase_ICA+GSR", name);
    % KS.phase_ICA_GSR_WM_CSF_MT_CEN = draw_FCD_pdf_KS(FCD_phase, "phase_ICA+GSR+WM_CSF_MT_CEN", name);

    % pack into metrics
    metrics = struct();
    metrics.FC_corr = FC_corr;
    metrics.KS = KS;
end
