function [W1,W2,wlparam,WCo,cl,S1,S2,S12]=wl_cohSingleTrial(dat1,dat2,dt,df,frange,opt,mother,wlopt,nco,ncy);
%function [W1,W2,wlparam,WCo,cl,[S1],[S2],[S12]]=wl_cohSingleTrial(dat1,dat2,dt,df,frange,opt,mother,wlopt,nco,ncy);
%
% Single trial coherence includes wavelet transforms of the two time series
%
% Parameters
%   x           Time series
%   dt          Time-resolution
%   df          Frequency-resolution (as a Fourier period)
%   frange      (opt) Frequency range (may be empty)
%   opt         Analysis options
%                   p - zero pad to power of 2
%                   l - Linear frequency scale (frange required)
%                   w - Fixed integration window (default: freq dependent)
%   mother      Mother wavelet as a string
%   wlopt       Options for the mother wavelet
%               (Passed as a cell array, i.e.: {1,2,3,...})
%   nco         No. oscillations in wavelet
%   ncy         No. oscillations in integration window
%
% NB: Morlet specific implementation (relating to w0 calculation)
%
%function [W1,W2,wlparam,WCo,cl,[S1],[S2],[S12]]=wl_cohSingleTrial(dat1,dat2,dt,df,frange,opt,mother,wlopt,nco,ncy);

% Default parameters
spectout=0;
crossout=0;
details=0;

% Check input parameters
if (nargin<10)
    error('Not enough input parameters');
end;

% Check output parameters
if (nargout>=6)
    spectout=1;
end;
if (nargout>=8)
    crossout=1;
end;

% Parse options string
options=deblank(opt);
opt_str='';
bicoh_opt='';
while (any(options))                % Parse individual options from string.
	[opt,options]=strtok(options);
	optarg=opt(2:length(opt));      % Determine option argument.
	switch (opt(1))
        case 'd'                    % Display progress
            details=1;
        case 'w'                    % Fixed integration window
            bicoh_opt=[bicoh_opt ' ' opt];
    	otherwise
            opt_str=[opt_str ' ' opt];
    end;
end;

% Require wavelet to produce n_co oscilations
w0=2*pi;
wlopt={w0};

% Transform the two signals
if (details), disp('Wavelet Transform Data Channel 1'); end;
[W1,wlparam]=wlfTransform(dat1,dt,df,frange,opt_str,mother,wlopt);
if (details), disp('Wavelet Transform Data Channel 2'); end;
[W2,wlparam]=wlfTransform(dat2,dt,df,frange,opt_str,mother,wlopt);

% Calculate single trial coherence
if (details), bicoh_opt=[bicoh_opt ' d']; end;
if (spectout)
    if (crossout)
        [WCo,cl,S1,S2,S12]=wlCohSingleTrial(W1,W2,wlparam,1/dt,nco,ncy,bicoh_opt);
    else
        [WCo,cl,S1,S2]=wlCohSingleTrial(W1,W2,wlparam,1/dt,nco,ncy,bicoh_opt);
    end;
else
    [WCo,cl]=wlCohSingleTrial(W1,W2,wlparam,1/dt,nco,ncy,bicoh_opt);
end;
