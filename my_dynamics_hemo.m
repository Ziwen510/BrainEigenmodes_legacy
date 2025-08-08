clearvars;         
close all force;
addpath(genpath('functions_matlab'));

function tc = load_tc(filename)
%LOAD_TC   Load the time-course matrix 'tc' from a MAT-file.
%   tc = LOAD_TC(filename)  
%   - First attempts S = load(filename,'tc'); tc = S.tc;
%   - If that errors (e.g. data >4GB in v7 MAT), it will try h5read(filename,'/tc').
%
%   Example:
%     tc = load_tc('my_mats/cortical_str_nuee0.2000_nues0.1000.mat');

    % validate input
    if ~isfile(filename)
        error('File not found: %s', filename);
    end

    % try standard MAT-file load
    try
        S = load(filename, 'tc');
        tc = S.tc';
        return
    catch ME
        warning('Standard load failed: %s\nFalling back to HDF5…', ME.message);
    end

    % fallback to HDF5 read
    try
        tc = h5read(filename, '/tc');
    catch H5E
        error('HDF5 read failed for %s:\n%s', filename, H5E.message);
    end
end

surface_interest = 'fsLR_32k';
hemisphere = 'lh';
mesh_interest = 'midthickness';

[vertices, faces] = read_vtk(sprintf('data/template_surfaces_volumes/%s_%s-%s.vtk', surface_interest, mesh_interest, hemisphere));
surface_midthickness.vertices = vertices';
surface_midthickness.faces = faces';

% Load cortex mask
cortex = dlmread(sprintf('data/template_surfaces_volumes/%s_cortex-%s_mask.txt', surface_interest, hemisphere));
cortex_ind = find(cortex);


%% Tmax HERE
tmax = 865 * 1000;  % in ms
tstep = 720; % in ms

tspan = [0, tmax];
T = 0:tstep:tmax;
disp("original length of specified T");
disp(size(T));

hemisphere = 'lh';
input = 'noise';
Tmax = tmax / 1000;
simulator = 'cortical';
% norm = '_dianorm';
% norm = '_globalnorm';
norm = '';
version = '_newacc';
% version = '';
rmax = '430';
variable = '_bold';
% boldnorm = '_perframe';
boldnorm = '';
cropping = '';

rate = 1000;
nu_es = 0.0001 * rate;
% nu_ee = 0.00021238 * rate;
% v0 = 0.00339808 * rate;
% theta = 0.015 * rate;
% sigma = 0.0033 * rate;
nu_ee = 0.00006 * rate;
v0 = 0.0006 * rate;
theta = 0.0126766 * rate;
sigma = 0.0038 * rate;
r_s = 18;
v = r_s * 116;

switch input
    case 'noise'
        seed = 9;
        mean_I = 0;
        ASD = 0.01;
        str = sprintf("white_ASD%.4f_seed%d_mean%.1f", ASD, seed, mean_I);
    case 'pulse'
        cite = 13289;
        % cite = 18245;
        amp = 1000;
        str = sprintf("pulse_%d_I%d", cite, amp);
end
S = load_tc(sprintf('my_mats/cortical%s_%s_nuee%.4f_nues%.4f_r%d_v%d_v0%.4f_theta%.4f_sigma%.4f_T%.0f_rmax%s%s%s%s.mat', version, str, nu_ee, nu_es, r_s, v, v0, theta, sigma, Tmax, rmax, norm, boldnorm, variable));   % loads only the variable 'tc' into the struct S
fname = sprintf("my_%s%s_%s_r%d_T%.1f%s%s", simulator, version, str, r_s, Tmax, norm, cropping);
disp(fname);

simulated_activity_rest = S;
disp("Size of loaded TC:")
disp(size(simulated_activity_rest));


warmup = 40000; % in ms
target_tstep = 720;
switch cropping
    case "_cropped"
        T_target = warmup:target_tstep:217*1000;                         % [1 × M]
    otherwise
        T_target = warmup:target_tstep:tmax;
end

idx = (round(T_target./tstep) + 1);                     % +1 because T(1)=0 at index 1
tc_down = simulated_activity_rest(:, idx);                                 % [Ncortex × M']


% simulated_actitvity_rest is [Ncortex × T]
[nc, Tn] = size(tc_down);
Ntotal    = numel(cortex);
full_activity = zeros(Ntotal, Tn);
full_activity(cortex_ind, :) = tc_down;

disp("size of final full tc:")
disp(size(full_activity));

metrics = analysis_pipeline(full_activity, fname, T_target, cortex, surface_midthickness);

disp(metrics);