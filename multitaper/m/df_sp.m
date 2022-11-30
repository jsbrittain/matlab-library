function [DF11,DF22,DF12,params]=df_sp(dat1,dat2,rate,opt_str,v,dk)
%function [DF11,DF22,DF12,params]=df_sp(dat1,dat2,rate,opt_str,[v,dk])
%
% Dual-frequency analysis
%
% Input parameters
%       dat1        Time series 1
%       dat2        Time series 2
%       rate        Sampling rate
%       opt_str     Options string
%                       r<0,1,2>    rectify (none, ch1, ch1&2)
%                       W           multitaper bandwidth (default: 6, 0=none)
%                       f           maximum analysis frequency (default: Nyquist)
%                       p           phase coherence
%                       a           amplitude coherence
%       w           (opt) Data tapers
%       dk          (opt) Associated taper weights
%
%function [DF11,DF22,DF12,params]=df_sp(dat1,dat2,rate,opt_str,[v,dk])

% Check input parameters
if ((nargin~=4) & (nargin~=6))
    error(' Incorrect number of input parameters');
end;

% Default parameters
W=6;                        % multitaper bandwidth
maxf=rate/2;                % maximum analysis frequency (Nyquist)
multitaper=true;            % perform multitaper analysis
phasecoherence = false;
amplitudecoherence = false;
phaseamplitude = false;

% Parse options string
options=deblank(opt_str);
opt_str='';
while (any(options))                % Parse individual options from string.
	[opt,options]=strtok(options);
	optarg=opt(2:length(opt));      % Determine option argument.
	switch (opt(1))
	    case 'r'                    % Rectify channels
        	n=str2num(optarg);
			if ((n<0) || (n>2))
			    error(['Error in option argument -- r' optarg]);
			end;  
			if (n~=1)   % Rectify ch 1.
			    dat1=abs(dat1-mean(dat1));
			end;  
			if (n>=1)   % Rectify ch 2.
			    dat2=abs(dat2-mean(dat2));
			end;
        case 'W'
            W=str2num(optarg);
            if (isempty(W))
                error(['Error in option argument -- W' optarg]);
            end;
            if (W==0)
                multitaper=0;
            end;
        case 'p'
            phasecoherence = true;
        case 'a'
            amplitudecoherence = true;
        case 'q'
            phaseamplitude = true;
        case 'f'
            maxf=str2num(optarg);
            if (isempty(W))
                error(['Error in option argument -- f' optarg]);
            end;
        otherwise                   % Options for wavelet analysis
            error(['Unknown option -- ' opt]);
    end;
end;

% Get data parameters
N=length(dat1);

% Choose tapers
if (~exist('v'))
	% Generate Discrete Prolate Spheroidal Sequences (DPSS)
	if (multitaper)
        [v,dk]=dpss(N,W);     % 3rd part routine (Breitenberger)
	else
        v=ones(N,1);
	end;
    % Determine weighting scheme
    dk=dk/sum(dk);
end;
% Determine no. tapers
tapercount=size(v,2);

% Determine frequency range and maximum frequency
freqs=rate*[0:N/2]/N;
fmax=min(floor(N/2),find(abs(freqs-maxf)==min(abs(freqs-maxf)))); % pts
freqs=freqs(1:fmax);

% Allocate memory
DF11=zeros(fmax);
DF22=zeros(fmax);
DF12=zeros(fmax);

% Zero mean data
dat1=dat1-mean(dat1);
dat2=dat2-mean(dat2);

% Form multitaper estimate
for taper=1:tapercount
    % Fourier Transform tapered data and truncate
    F1=fft(v(:,taper).*dat1); F1=F1(1:fmax);
    F2=fft(v(:,taper).*dat2); F2=F2(1:fmax);
    % Normalise if phase coherence
    if ( phasecoherence )
        F1 = F1./abs(F1);
        F2 = F2./abs(F2);
    end;
    if ( amplitudecoherence )
        F1 = abs(F1);
        F2 = abs(F2);
    end;
    % Calculate dual-frequency spectra
    DF11=DF11+dk(taper)*(F1*F1');               % NB: Conjugate transpose
    DF22=DF22+dk(taper)*(F2*F2');
    DF12=DF12+dk(taper)*(F2*F1');
end;

% Determine equivalent no. segments
L=1/sum(dk.^2);

% Assign output parameters
params.W=W;
params.freqs=freqs;
params.fmax=maxf;
params.tapercount=tapercount;
params.L=L;
