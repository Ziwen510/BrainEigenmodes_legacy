function BOLD = model_BOLD_balloon_vertex(ext_input, param, method)
% model_BOLD_balloon_vertex.m
% 
% Simulate the BOLD balloon model directly on the cortical surface (vertex level)
% without using eigenmodes.
% 
% Inputs:
%   ext_input : [V x T] external input time series per vertex
%   param     : struct of model parameters
%               Required fields:
%                 T     : [1 x T] time vector
%                 tstep : scalar time step
%                 V0    : resting blood volume fraction
%                 k1,k2,k3 : BOLD coefficients
%                 kappa, gamma, tau, alpha, rho : hemodynamic parameters
%   method    : string, 'ODE' (default) or 'Fourier'
% 
% Outputs:
%   state     : struct containing fields z,f,v,q each [V x T]
%   BOLD      : [V x T] simulated BOLD signal per vertex

if nargin<3 || isempty(method)
    method = 'ODE';
end

[V, T] = size(ext_input);

switch upper(method)
    case 'ODE'
        % initialize the 4‐state F = [z f v q] at t=1 (double)
        F = repmat([0, 1, 1, 1], V, 1);

        % preallocate full BOLD in double
        BOLD = zeros(V, T);
        BOLD(:,1) = 100*param.V0 * ( ...
            param.k1*(1 - F(:,4)) + ...
            param.k2*(1 - F(:,4)./F(:,3)) + ...
            param.k3*(1 - F(:,3)) );

        % time‐stepping: never store z,f,v,q histories
        for k = 2:T
            S  = ext_input(:, k-1);              % [V×1]
            dF = balloon_ODE(S, F, param);       % returns [V×4] double
            F  = F + param.tstep * dF;           % update current state

            % compute and store BOLD(:,k)
            BOLD(:,k) = 100*param.V0 * ( ...
                param.k1*(1 - F(:,4)) + ...
                param.k2*(1 - F(:,4)./F(:,3)) + ...
                param.k3*(1 - F(:,3)) );

            % clear temporaries (optional)
            % clear S dF
        end
        
    case 'FOURIER'
        %--- Fourier-based simulation at vertex level using balloon_Fourier ---
        % 1) build symmetric time vector around zero
        T_append = -param.T(end):param.tstep:param.T(end);
        Nt       = numel(T_append);
        % find index corresponding to t = 0
        t0_ind   = dsearchn(T_append', 0);

        % 2) pad each vertex’s input with zeros for negative times
        ext_input_pad = zeros(V, Nt);
        ext_input_pad(:, t0_ind:end) = ext_input;

        % 3) preallocate the full padded output
        BOLD_pad = zeros(V, Nt);

        % 4) for each vertex, call the same balloon_Fourier routine
        for v = 1:V
            coeff = ext_input_pad(v, :);           % [1×Nt]
            yout  = balloon_Fourier(coeff, T_append, param);
            BOLD_pad(v, :) = yout;                 % [1×Nt]
        end

        % 5) extract only the original (non-negative-pad) segment
        BOLD  = BOLD_pad(:, t0_ind:end);
        
    otherwise
        error('Unknown method "%s". Choose "ODE" or "Fourier".', method);
end
end


function dF = balloon_ODE(S, F, param)
% balloon_ODE.m
%
% Calculate the temporal activity by solving an ODE.
%
% Inputs: S          : spatiotemporal external input [V x 1]
%         F          : solutions at one time point
%         param      : model parameters (struct)
%
% Output: dF         : time derivative of variables [4 x 1]
%
% Original: James Pang, Monash University, 2022

z = F(:,1);
f = F(:,2);
v = F(:,3);
q = F(:,4);

dF(:,1) = S - param.kappa*z - param.gamma*(f - 1);
dF(:,2) = z;
dF(:,3) = (1/param.tau)*(f - v.^(1/param.alpha));
dF(:,4) = (1/param.tau)*((f/param.rho).*(1 - (1 - param.rho).^(1./f)) - q.*v.^(1/param.alpha - 1));

end

function out = balloon_Fourier(mode_coeff, T, param)
% balloon_Fourier.m
%
% Calculate the temporal activity of one mode via Fourier transform.
%
% Inputs: mode_coeff : coefficient of the mode [1 x T]
%         T          : time vector with zero center [1 x T]
%         param      : model parameters (struct)
%
% Output: out        : activity [1 x T]
%
% Original: James Pang, Monash University, 2022

Nt = length(T);
Nw = Nt;
wsamp = 1/mean(param.tstep)*2*pi;
jvec = 0:Nw-1;
w = (wsamp)*1/Nw*(jvec - Nw/2);

% mode_coeff_fft = ctfft(mode_coeff, param.T);	
mode_coeff_fft = coord2freq_1D(mode_coeff, w);
	
% Transfer functions from Robinson et al. 2006, Aquino et al. 2012, 2014,
% Pang et al. 2016, 2018
T_Fz = 1 ./ (-(w + 1i*0.5*param.kappa).^2 + param.w_f^2);
T_yF = param.V_0 * (param.alpha*(param.k2 + param.k3)*(1 - 1i*param.tau*w) - (param.k1 + param.k2)*(param.alpha + param.beta - 1 - 1i*param.tau*param.alpha*param.beta*w))./((1 - 1i*param.tau*w).*(1 - 1i*param.tau*param.alpha*w));
T_yz = T_yF.*T_Fz;

out_fft = T_yz.*mode_coeff_fft;

% calculate inverse Fourier transform
out = real(freq2coord_1D(out_fft, w));

end