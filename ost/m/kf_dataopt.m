function [dat1,dat2,opt_str]=kf_dataopt(dat1,dat2,opt_str);
%function [dat1,dat2,opt_str]=kf_dataopt(dat1,dat2,opt_str);
% Parse options string and return processed data
%
% Options
%       a<len>      High pass filter data
%       r<0|1|2>    Rectify channel 1, 2 or both
%

% Default options
filt_no=1;
filt(1).windowSize=0;
filt(2).windowSize=0;
rectify(1)=logical(0);
rectify(2)=logical(0);

% Parse options string
options=deblank(opt_str);
opt_str=[];
while (any(options))                            % Parse individual options from string.
	[opt,options]=strtok(options);
	optarg=opt(2:length(opt));                  % Determine option argument.
	switch (opt(1))
	    case 'a'                                % Channel reference invert option.
            filt(filt_no).windowSize=str2num(optarg);
            filt_no=2;
        case 'r'
            switch (str2num(optarg))
                case 0, rectify(1)=1;
                case 1, rectify(2)=1;
                case 2
                    rectify(1)=1;
                    rectify(2)=1;
            end;
	    otherwise
	        opt_str=[opt_str ' ' opt(1) optarg];
	end
end

% Return any remaining options
if (isempty(opt_str))
    opt_str=' ';
end;

% Perform low-pass filtering
if (filt(1).windowSize~=0)
    dat1=dat1-filter(ones(1,filt(1).windowSize)/filt(1).windowSize,1,dat1);
end;
if (filt(2).windowSize~=0)
    dat2=dat2-filter(ones(1,filt(2).windowSize)/filt(2).windowSize,1,dat2);
end;

% Rectify channels
if (rectify(1))
    dat1=abs(dat1);
end;
if (rectify(2))
    dat2=abs(dat2);
end;
