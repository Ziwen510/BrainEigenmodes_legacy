function pipeline_from_cortexBOLD(TC, cortex_id, fname, start_point)
    V_full = 32492;

    if ~isnumeric(TC)
        mfile = matfile(TC,'Writable',false);
        [V, T_full] = size(mfile, 'simulated_activity_rest');
        fprintf('Full TC dims: %d vertices × %d timepoints\n', V, T_full);
    
        % Preallocate only the cortex‐subset
        % nC = numel(cortex_id);
        % sample = mfile.simulated_activity_rest(1,1);
        % TC     = zeros(nC, T_full, 'like', sample);
        % TC = mfile.simulated_activity_rest(cortex_id, :);
    
        TC = mfile.simulated_activity_rest;
        clear mfile
        disp('Loaded TC');
    end
    
    %% 7) Crop warmup
    TC_balloon = TC(:, start_point:end);
    croppedT   = size(TC_balloon,2);
    disp('Cropping done. Size after cropping:');
    disp(size(TC_balloon));
    
    %% 8) Re-embed into full space only when needed
    TC_full = zeros(V_full, croppedT, 'like', TC_balloon);
    TC_full(cortex_id, :) = TC_balloon;
    clear TC_balloon
    disp('Full TC done');

    %% 9) Compute FC and save
    hemisphere = 'lh';
    data_parc = my_parcellate(TC_full, hemisphere);


    % FC_vec = draw_FC(data_parc, sprintf('%s_FC', fname));
    % FCD = draw_FCD(data_parc, 83, sprintf('%s_FCD', fname));
    draw_FCD_phase(data_parc, fname);
    clear TC_full
    disp('Pipeline complete');
end
