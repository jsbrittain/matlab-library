function [dat1,dat2,opt_str]=kf_dataopt(dat1,dat2,opt_str);
%function [dat1,dat2,opt_str]=kf_dataopt(dat1,dat2,opt_str);
% Parse options string and return processed data
%
% Currently only performs low pass filtering (option 'a')
%

% Default options
filt_no=1;
filt(1).windowSize=0;
filt(2).windowSize=0;

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
