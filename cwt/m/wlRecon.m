function x=wlRecon(W,wlparam,srange);
%function x=wlRecon(W,wlparam,[srange]);
%
% Reconstruct a time series from a non-orthogonal basis
% Suitable for use with real-valued time series only
% (For applicability to complex time-series remove the
%  real() constrain in the main equation)
%
% Input parameters
%   W       Wavelet transform
%   wlparam Wavelet transform parameters
%   srange  (opt) Scale range for reconstruction
%
% Output
%   x       Time series
%

% Determine input parameters
bandpass=0;
if (nargin>2)   % If scale range provided
    if (~isempty(srange))   % If srange is not empty
        if (wlparam.linearscale)    % If a linear scale
            scale=srange;
        else                        % If a base 2 scale
            s0=srange(1);
            J=srange(2);
        end;
    else                    % If srange empty
        if (wlparam.linearscale)    % If a linear scale
            scale=wlparam.scale';
        end;
    end;
    bandpass=1;
end;

% Obtain data parameters
N=size(W,2);

% Determine wavelet parameters
mother=wlparam.mother;
wlopt=wlparam.wlopt;
if (~wlparam.linearscale)
    dj=wlparam.dj;
else
    dj=wlparam.df/wlparam.fourier_factor;   % BODGE - FIX THIS
end;
dt=wlparam.dt;

% Determine scale (If linearscale then its already defined)
if (~wlparam.linearscale)
	if (~bandpass)
		% Calculate scale range
		s0=2*dt;
		J=ceil((1/dj)*log2(N*dt/s0));
	end;
	j=0:J;
	scale=s0.*2.^(j*dj)';
else
    if (~bandpass)
        scale=wlparam.scale';
    end;
end;

% Reduce W to that of 'scale'
if (bandpass)
    if (wlparam.linearscale)
        % Calculate full scale as relating to W (linear scale)
        Wscale=wlparam.scale';
    else
        % Calculate full scale as relating to W (base 2 scale)
        Ws0=2*dt;
		WJ=ceil((1/dj)*log2(N*dt/Ws0));
        Wj=0:WJ;
        Wscale=Ws0.*2.^(Wj*dj)';
    end;
    
    % Reduce W scale size from Wscale to that of 'scale'
    WSstart=dsearchn(Wscale,scale(1));  % Closest point search
    WSend=dsearchn(Wscale,scale(length(scale)));
    W=W(WSstart:WSend,:);
end;

% Determine frequency range
halflength=floor(N/2);
omega=[(2*pi*(0:halflength)/(N*dt)) (-(2*pi*((halflength+1):(N-1)))/(N*dt)) ]';

% Obtain wavelet parameters
wl_mother=['wl_' mother];
if (~exist('wlparam.Cg'))
    wlp_mother=['wlp_' mother];
    tmp_param=feval(wlp_mother,wlparam.omega,wlparam.scale,wlparam.dt,dj,wlparam.wlopt{:});
    wlparam.Cg=tmp_param.Cg;
end;
psi0=feval(wl_mother,0,0);

% Reconstruct time-series
for n=1:N
    x(n,1)=(dj*sqrt(dt)/(wlparam.Cg*psi0))*sum(real(W(:,n))./sqrt(scale));
end;
