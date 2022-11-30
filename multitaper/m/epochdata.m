function [epochs,erptime]=epochdata(dat,trig,rate,offset,duration,baseline,opt_str)
%function [epochs,erptime]=epochdata(dat,trig,rate,offset,duration,baseline,opt_str)
%
% Construct epoch data
%
% Options string
%       b<n>        0: subtract mean baseline per trial
%                   1: subtract trimmean (over trials) of baseline means (over time)
%                   2: none
%
%function [epochs,erptime]=epochdata(dat,trig,rate,offset,duration,baseline,opt_str)

% Check input parameters
if (~exist('opt_str','var'))
    opt_str = '';
end;

% Define default parameters
baseline_mode = 0;  % 0: individual, 1: mean, 2=none

% Parse options string
options=deblank(opt_str);
while (any(options))                % Parse individual options from string.
	[opt,options]=strtok(options);
	optarg=opt(2:length(opt));      % Determine option argument.
	switch (opt(1))
	    case 'b'        % Baseline mode
            baseline_mode = str2num(optarg);
            if (isempty(baseline_mode))
                error(['Error in option argument -- b' optarg]);
            end;
            if ((baseline_mode<0) || (baseline_mode>2))
                error(['Invalid option in argument -- b' optarg]);
            end;
        otherwise                   % Unknown options
            error(['Unknown option -- ' opt]);
    end;
end;

% Determine epoch timespan
tt = round(offset*rate/1000)+[0:duration*rate/1000];     % samples

% Remove trials exceeding epoch limits
trig = trig(((trig+tt(1))>0) & ((trig+tt(end))<=length(dat)) & (~isnan(trig)));

% Reserve memory
epochs = zeros(length(tt),size(dat,2),length(trig));

% Construct epochs
for ind = (1:length(trig))
    epochs(:,:,ind)=dat(trig(ind)+tt,:);
end;

% Output time vector
erptime = tt*1000/rate;  % msecs

% Baseline
if (exist('baseline','var'))
    if (~isempty(baseline))
        % Get baseline region
        tspan=dsearchn(erptime',baseline'); tspan=(tspan(1):tspan(2));
        % Recurse trials and remove baseline
        for ch=(1:size(epochs,2))
            % Determine baseline correction mode
            switch ( baseline_mode )
                case 2,    % None
                case 0,     % Baseline per trial
                    for ind=(1:size(epochs,3))
                        epochs(:,ch,ind) = epochs(:,ch,ind) - mean(epochs(tspan,ch,ind),1);
                    end;
                case 1,     % Determine mean baseline over time & trials then subtract from each epoch
                    blcorrection = trimmean( mean(epochs(tspan,ch,:),1), 5, 3);     % Trim mean
                    epochs(:,ch,:) = epochs(:,ch,:) - blcorrection;
            end;
        end;
    end;
end;

% Determine output parameters
if (nargout<2)
    clear('erptime');
end;
