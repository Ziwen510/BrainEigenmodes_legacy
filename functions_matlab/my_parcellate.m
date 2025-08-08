function data_parc = my_parcellate(data, hemisphere)
    parc_name = 'Glasser360';
    parc = dlmread( sprintf('data/parcellations/fsLR_32k_%s-%s.txt', ...
                            parc_name, hemisphere) );
    
    %--- parcellate & normalize timeseries -------------------------
    data_parc = calc_parcellate(parc, data);        % [parcel × time]
    % data_parc = calc_normalize_timeseries(data_parc');  % [time × parcel]
    data_parc(isnan(data_parc)) = 0;
end