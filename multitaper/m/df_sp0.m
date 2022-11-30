function [df11,df22,df12,params,coh11,coh22,coh12]=df_sp0(dat1,dat2,duration,rate,opt_str)
%[df11,df22,df12,params,[coh11,coh22,coh12]]=df_sp0(dat1,dat2,duration,rate,opt_str)
%
% Multitaper analysis
% Type 0 by default
%
%
% Options as mt_sp2
%
%[df11,df22,df12,params,[coh11,coh22,coh12]]=df_sp0(dat1,dat2,duration,rate,opt_str)

% Default parameters
dt = [];

% Parse options string
options=deblank(opt_str);
opt_str='';
while (any(options))                % Parse individual options from string.
	[opt,options]=strtok(options);
	optarg=opt(2:length(opt));      % Determine option argument.
	switch (opt(1))
        case 't',                   % Time step for overlapping segments
            dt = str2num(optarg);
            if (isempty(dt))
                error(['Error in option argument -- W' optarg]);
            end;
            dt = round( dt*rate/1000 );
        otherwise                   % Options for segmented analysis
            opt_str=[opt_str ' ' opt];
    end;
end;

% Determine segments
offset=0;                       % Default offset
if (isempty(dt))
    dt=duration*rate/1000;          % msecs -> samples
end;
trig=[1:dt:length(dat1)-ceil(duration*rate/1000)];    % Disjoint sections

% Perform analysis over segments
[df11,df22,df12,params,coh11,coh22,coh12]=df_sp2(dat1,dat2,trig,offset,duration,rate,opt_str);
