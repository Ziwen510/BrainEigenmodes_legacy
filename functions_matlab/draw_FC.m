function FCvec = draw_FC(data_parc, fname)
    c_map = load("Rapaeh_color_table_FCD.mat", "c3");
    c_map = c_map.c3;

    data_parc = data_parc';
    % input: N*T
    FC_emp = corrcoef(data_parc);

    
    %--- plot ------------------------------------------------------
    hFig = figure;  
    imagesc(FC_emp);
    axis square
    colorbar
    clim([-1 1]); 
    colormap(c_map);
    title(fname,'Interpreter','none')
    xlabel('Parcel'), ylabel('Parcel')
    
    %--- ensure output folder exists & save ------------------------
    savePath = sprintf("./outputs/%s.png", fname);
    [outDir,~,~] = fileparts(savePath);
    if ~isempty(outDir) && ~exist(outDir,'dir')
        mkdir(outDir)
    end
    print(hFig, savePath, '-dpng', '-r300');
    
    %--- get upper‐triangle indices -------------------------------
    num_parcels = size(data_parc, 2);
    triu_ind = calc_triu_ind( zeros(num_parcels, num_parcels) );
    %--- return upper‐triangle vector -----------------------------
    FCvec = FC_emp(triu_ind);
    matSavePath = sprintf("./outputs/%s_FC.mat", fname);
    save(matSavePath, 'FCvec');
end
