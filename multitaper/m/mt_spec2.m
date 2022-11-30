function [S11,S22,S12,params]=mt_spec2(dat1,dat2,trig,offset,duration,rate,segdur,opt_str);
%function [S11,S22,S12,params]=mt_spec2(dat1,dat2,trig,offset,duration,rate,segdur,opt_str);
%
% Multi-taper spectrogram analysis type 2
%
% Input parameters
%   Data parameters
%       dat1        Time series 1
%       dat2        Time series 2
%       trig        Trigger
%       offset      Offset (ms)
%       duration    Duration (ms)
%       rate        Sampling rate
%       segdur      Segment duration
%   Analysis parameters
%       frange      Frequency range (may be [] for defaults)
%       opt_str     Options string (as mt_spectrogram.m)
%
%function [S11,S22,S12,params]=mt_spec2(dat1,dat2,trig,offset,duration,rate,segdur,opt_str);

% Determine data parameters
offset=round(offset*rate/1000);        % Convert ms -> samples
duration=round(duration*rate/1000);    %  |
segN=round(segdur*rate/1000);          %  |
trig=round(trig);
trig=trig(((trig+offset)>0) & ((trig+offset+duration)<length(dat1)));
trialcount=length(trig);

% Default parameters
NW=3;                   % time-bandwidth product

% Parse options string
options=deblank(opt_str);
while (any(options))                % Parse individual options from string.
	[opt,options]=strtok(options);
	optarg=opt(2:length(opt));      % Determine option argument.
	switch (opt(1))
        case 'W'        % Specify bandwidth
            NW=str2num(optarg);
            if (isempty(NW))
                error(['Error in option argument -- W' optarg]);
            end;
    end;
end;

% Generate DPSS sequences
[v,dk]=dpss(segN,NW);

% Recurively form trial-average spectrograms
for trial=1:trialcount
    % Display progress
    disp(['Multi-taper spectrogram trial ' int2str(trial) ' of ' int2str(trialcount)]);
    % Perform spectrogram anlaysis on single block
    trange=trig(trial)+offset+[0:duration-1];
    [S11a,S22a,S12a,params]=mt_spectrogram(dat1(trange),dat2(trange),rate,segdur,opt_str,v,dk);
    % Recursively form spectrograms
    if (trial==1)
        S11=S11a/trialcount;
        S22=S22a/trialcount;
        S12=S12a/trialcount;
    else
        S11=S11+S11a/trialcount;
        S22=S22+S22a/trialcount;
        S12=S12+S12a/trialcount;
    end;
end;

% Modify parameters structure
params.L=params.L*trialcount;
params.offset=offset*1000/rate;     % msecs
params.duration=duration*1000/rate; % msecs
