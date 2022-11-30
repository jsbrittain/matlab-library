function [sp,freqs] = phase_consistency( dat, rate, window, freqs, opt_str )
%
% Complex phase spectrum for a single channel
%
% FFT method
%

%error(' Does not work! Resets phase at each window!');
error(' Works, but only one sample per window so a bad approach!');

% Check input parameters
if (~exist('opt_str'))
    opt_str = '';
end;
if (isempty(window))
    window = length(dat)/rate;
end;

% Parameter defaults
maxf = rate/2;          % Maximum analysis frequency
display = true;

% Reshape input
dat = dat(:);
wnd = round(window*rate);       % samples
trials = floor(length(dat)/wnd);
dat = dat(1:trials*wnd);

% Parse options string
options=deblank(opt_str); opt_str='';
while (any(options))                % Parse individual options from string.
	[opt,options]=strtok(options);
	optarg=opt(2:length(opt));      % Determine option argument.
	switch (opt(1))
        case 'x'        % Display off
            display = false;
        otherwise
            error(['Unknown option -- ' opt]);
    end;
end;

% Recurse trials
sp = zeros(length(freqs),trials);
for k = (1:length(freqs))
    if (display)
        disp([' Evaluating frequency ' num2str(k) '/' num2str(length(freqs))]);
    end;
    for n = (1:trials)
        tt = (n-1)*wnd+(1:wnd).';
        % use tt for complex sinusoid to ensure wave always starts at t=0 so is consistent across time
        sp(k,n) = sum( (dat(tt)-mean(dat(tt))).*exp(1i*(2*pi*freqs(k)).*tt/rate) );
    end;
end;
