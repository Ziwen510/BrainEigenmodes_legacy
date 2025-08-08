function test_pipeline(TC, cortex_id, fname, tmax)
    % tc_matpath : string, path to .mat file containing variable full_tc
    % cortex_id  : vector of vertex-indices to keep
    % fname      : base name for output files
    % tmax       : simulation length in ms

    V_full = 32492;
    %% 1) Set up hemodynamic parameters
    bal_param = loadParameters_balloon_func;
    bal_param.tstep = 18 * 1e-3;
    bal_param.tmax  = tmax * 1e-3;
    bal_param.tspan = [0 bal_param.tmax];
    bal_param.T     = 0:bal_param.tstep:bal_param.tmax;

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

    

    %% 3) In-place frame-wise demeaning
    mu = mean(TC,1);        % 1×T
    % mu = 10;
    disp(mu);
    TC = TC - mu;           
    clear mu
    disp('Demeaning done');
    % disp(bal_param.T);
    

    %% 4) In-place global scaling
    % two-stage max avoids creating abs(TC(:))
    % colmax = max(abs(TC),[],1);
    % m      = max(colmax);
    % disp(m);
    % % scale  = 0.01 / m;
    % % disp(scale);
    % % TC     = TC * scale;
    % % clear colmax m
    % disp('Scaling done');

    %% 6) Run the balloon model
    TC_balloon = model_BOLD_balloon_vertex(TC, bal_param, 'ODE');
    clear TC
    disp('Balloon ODE done');
    disp(size(TC_balloon));
    disp(size(bal_param.T));
    draw_tc_samples(TC_balloon(:,1:end), cortex_id, fname, bal_param.T(1:end), 10, false);
    

    %% 7) Crop warmup
    TC_balloon = TC_balloon(:, 301:end);
    croppedT   = size(TC_balloon,2);
    disp(size(TC_balloon));
    disp('Cropping done');
    

    %% 8) Re-embed into full space only when needed
    TC_full = zeros(V_full, croppedT, 'like', TC_balloon);
    TC_full(cortex_id, :) = TC_balloon;
    disp(size(TC_full));
    clear TC_balloon
    disp('Full TC done');

    %% 9) Compute FC and save
    hemisphere = 'lh';
    FC_vec = draw_FC(TC_full, sprintf('%s_pipeline_FC', fname), hemisphere);
    clear TC_full
    disp('Pipeline complete');
end
