function [DF11,DF22,DF12,params,coh11,coh22,coh12]=df_sp2(dat1,dat2,trig,offset,duration,rate,opt_str)
%function [df11,df22,df12,params,[coh11,coh22,coh12]]=df_sp2(dat1,dat2,trig,offset,duration,rate,opt_str)
%
% Dual-frequency analysis type 2
%
% Input parameters
%   Data parameters
%       dat1        Time series 1
%       dat2        Time series 2
%       trig        Trigger
%       offset      Offset (ms)
%       duration    Duration (ms)
%       rate        Sampling rate
%   Analysis parameters
%       opt_str     Options string
%                       r<0,1,2>    rectify (none, ch1, ch1&2)
%                       W           multitaper bandwidth (default: 6, 0=none)
%                       f           maximum analysis frequency (default: Nyquist)
%                       c<p,a,c,q>  use correlation measure
%                                      p - phase coherence         [ OK ]
%                                      a - amplitude coherence     [ WHAT DOES THIS MEAN? ]
%                                      c - amplitude correlation   [ CORRELATION OVER WINDOW LENGTH ]
%                                      q - phase-amplitude coupling
%                                    coh11, coh22, coh12 contain the amplitude correlations
%                                    instead of coherence as derived from the multitaper periodograms
%
%function [df11,df22,df12,params,[coh11,coh22,coh12]]=df_sp2(dat1,dat2,trig,offset,duration,rate,opt_str)

% Default parameters
W=6;                        % multitaper bandwidth
maxf=rate/2;                % maximum analysis frequency (Nyquist)
multitaper=true;            % perform multitaper analysis
correlation=false;

% Parse options string
options=deblank(opt_str);
opt_str='';
while (any(options))                % Parse individual options from string.
	[opt,options]=strtok(options);
	optarg=opt(2:length(opt));      % Determine option argument.
	switch (opt(1))
        case 'W'
            W=str2num(optarg);
            if (isempty(W))
                error(['Error in option argument -- W' optarg]);
            end;
            if (W==0)
                multitaper=0;
            end;
        case 'f'
            maxf=str2num(optarg);
            if (isempty(maxf))
                error(['Error in option argument -- f' optarg]);
            end;
            opt_str=[opt_str ' ' opt];
        case 'c'                    % Correlation measures
            switch (optarg)
                case 'p', opt_str=[opt_str ' p'];
                case 'a', opt_str=[opt_str ' a'];
                case 'c', correlation = true; amporphase = [ 1 1 ];
                case 'q', correlation = true; amporphase = [ 0 1 ];
                otherwise
                    error(' Unknown correlation measure requested!');
            end;
        otherwise                   % Options for segmented analysis
            opt_str=[opt_str ' ' opt];
    end;
end;

% Determine data parameters
offset=offset*rate/1000;        % Convert ms -> samples
duration=duration*rate/1000;    %
trig=round(trig);
trialcount=length(trig);

% Generate Discrete Prolate Spheroidal Sequences (DPSS)
if (multitaper)
    if (W==0.5)
        [v,dk]=dpss(duration,1);     % 3rd part routine (Breitenberger)
        v=v(:,1); dk=dk(1);
    else
        [v,dk]=dpss(duration,W);     % 3rd part routine (Breitenberger)
    end;
else
    v=ones(duration,1); dk=1;
end;
tapercount=size(v,2);

% Determine weighting scheme
dk=dk/sum(dk);

% Determine frequency range and maximum frequency
freqs=rate*[0:duration/2]/duration;
fmax=min(floor(duration/2),find(abs(freqs-maxf)==min(abs(freqs-maxf)))); % pts
%fmax=dsearchn(freqs.',maxf);
freqs=freqs(1:fmax);

% Allocate memory
if ( ~correlation )
    DF11=zeros(fmax);
    DF22=zeros(fmax);
    DF12=zeros(fmax);
else
    DF11=zeros(fmax,trialcount);
    DF22=zeros(fmax,trialcount);
    DF12=zeros(fmax,trialcount);
end;

% Recursively form dual-frequency spectrum
for trial = (1:trialcount)
    % Determine time range for trial
    disp(['Dual-frequency trial ' int2str(trial) ' of ' int2str(trialcount)]);
    trange=trig(trial)+offset+(0:duration-1);
    % Perform single segment analysis
    if (~correlation)
        % Dual-frequency spectra
        [df11,df22,df12,params]=df_sp(dat1(trange),dat2(trange),rate,opt_str,v,dk);
        % Update average
        DF11=DF11+df11/trialcount;
        DF22=DF22+df22/trialcount;
        DF12=DF12+df12/trialcount;
    else
        % Complex spectra
        [df11,df22,df12,params]=mt_sp(dat1(trange),dat2(trange),rate,[opt_str ' c'],v,dk);
        if (amporphase(1))
            df11 = abs(df11);
        else
            df11 = angle(df11);
        end;
        if (amporphase(2))
            df22 = abs(df22);
        else
            df22 = angle(df22);
        end;
        % Standard spectra for correlations
        DF11(2:end,trial) = df11;
        DF22(2:end,trial) = df22;
        DF12(2:end,trial) = df12;
    end;
end;
% Amplitude correlation over window length
if ( correlation )
    disp('Calculating correlations...');
    %fn = inline('log10(abs(x))');
    fn = inline('x');
    coh11 = zeros(fmax);
    coh22 = zeros(fmax);
    coh12 = zeros(fmax);
    for n1 = (1:size(DF11,1))
        for n2 = (1:size(DF11,1))
            R = corrcoef( fn(squeeze(DF11(n1,:))), fn(squeeze(DF11(n2,:))) );
            coh11(n2,n1) = R(1,2);
            R = corrcoef( fn(squeeze(DF22(n1,:))), fn(squeeze(DF22(n2,:))) );
            coh22(n2,n1) = R(1,2);
            R = corrcoef( fn(squeeze(DF11(n1,:))), fn(squeeze(DF22(n2,:))) );
            coh12(n2,n1) = R(1,2);
        end;
    end;
    DF11 = mean(DF11,2);
    DF22 = mean(DF22,2);
    DF12 = mean(DF12,2);
end;

% Modify parameters structure
params.W=W;
params.freqs=freqs;
params.fmax=maxf;
params.tapercount=tapercount;
params.trialcount=trialcount;
params.L=trialcount*params.L;
params.correlation=correlation;

% Check output arguments
if ( (nargout>4) && (~correlation) )
    [coh11,coh22,coh12] = df_sp2_coh(DF11,DF22,DF12);
end;
