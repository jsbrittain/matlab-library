function x=wlfRecon(W,wlparam,frange);
%function x=wlfRecon(W,wlparam,[frange]);
%
% Reconstruct a time series from a non-orthogonal basis
% Suitable for use with real-valued time series only
% (For applicability to complex time-series remove the
%  real() constrain in the main equation)
%
% Input parameters
%   W       Wavelet transform
%   wlparam Wavelet transform parameters
%   frange  (opt) Frequency range for reconstruction
%
% Output
%   x       Time series
%

% Determine input parameters
bandpass=0;
if (nargin>2)
    if (~isempty(frange))
        f_min=frange(1);
        f_max=frange(2);
        bandpass=1;
    else
        f_min=df;
        f_max=1/(2*wlparam.dt);
    end;
end;

% Convert Fourier parameters to scale
dj=wlparam.df/wlparam.fourier_factor;

% Call reconstruction function
if (bandpass)
    if (wlparam.linearscale)
        % Linear scale
        srange=1./(f_min:wlparam.df:f_max);             % Inverse frequency range = period range
        srange=(srange/wlparam.fourier_factor)';           % Fourier period -> scale
    else
        % Frequency to scale conversion
        s_max=1./(wlparam.fourier_factor*f_min);
        s_min=1./(wlparam.fourier_factor*f_max);
        
        % Determine scale range
        s0=s_min;
        N=size(W,2);
        J=ceil((1/dj)*log2(s_max/(wlparam.fourier_factor*s0)));
        srange=[s0 J];
    end;
    
    % Reconstruct in terms of wavelet scale
    x=wlRecon(W,wlparam,srange);
else
    x=wlRecon(W,wlparam);
end;
