function wlScaloSummary(time,dat,power,wlparam,opt);
%function wlScaloSummary(time,dat,power,wlparam,opt);
%
% Plot raw data, scalogram & global wavelet power
% Returns having set the current axis to the center scalogram
%
% Input parameters
%   time        Time vector
%   dat         Time series
%   power       Transform power
%   wlparam     Wavlet transform parameters
%   opt         (opt) Options string
%                   l<n>    Plot power log10 (0=none, 1=scalogram, 2=scalogram and global power)
%                   i       Image plot (default: contour plot)
%
%function wlScaloSummary(time,dat,power,wlparam,opt);

% Routine written to support frequency plotting, therefore wavelet
% transforms performed using wlTransform have frequency content added
% determined before plotting as frequencies.

% Default parameters
logplot=1;
plotstyle=0;        % 0-contour, 1-image

% Check if options string is provided
if (nargin<5)
    opt='';
end;

% Parse options string
options=deblank(opt);
opt_str='';
while (any(options))                % Parse individual options from string.
	[opt,options]=strtok(options);
	optarg=opt(2:length(opt));      % Determine option argument.
	switch (opt(1))
        case 'l'                    % Log plotting
            logplot=str2num(optarg);
        case 'i'                    % Image plot
            plotstyle=1;
            opt_str=[opt_str ' ' opt];
    	otherwise
            opt_str=[opt_str ' ' opt];
    end;
end;

% Determine frequency parameters if original analysis was conducted in scale
if (~isfield(wlparam,'freqs'))
    wlparam.freqs=1./(wlparam.scale*wlparam.fourier_factor);
    wlparam.coi=[Inf (1./(wlparam.fourier_factor*wlparam.coi(wlparam.coi>0))) Inf];
end;

% Convert power to log10 scale (Traditional for spectra)
if (logplot==2)
    power=log10(power);
end;

% Calculate global wavelet spectrum
glPower=mean(power,2);

% Calculate global wavelet spectrum without COI region
if (wlparam.display_coi)
	power2=power;
	for ind=1:length(time)
        power2(wlparam.freqs<wlparam.coi(ind),ind)=0;
	end;
	glPower2=mean(power2,2);
end;

% Convert power to log10 scale (for spectrogram but not for global power)
if (logplot==1)
    power=log10(power);
end;

% Setup axes
figure;
ax(1)=subplot(3,1,1);
ax(2)=subplot(3,1,2);
ax(3)=subplot(3,1,3);

% Plot raw data
subplot(ax(1));
set(ax(1),'position',[.1 .85 .7 .1]);
plot(time,dat,'k-');
xlim([min(time) max(time)]);

% Scalogram
subplot(ax(2));
set(ax(2),'position',[.1 .1 .7 .7]);
scalogram(time,power,wlparam,opt_str);

% Global Wavelet power with and without COI
freqs=wlparam.freqs;
freqlim=[min(wlparam.freqs) max(wlparam.freqs)];
if (~wlparam.linearscale)
    freqs=log2(freqs);
    freqlim=log2(freqlim);
end;
subplot(ax(3));
set(ax(3),'position',[.85 .1 .1 .7]);
if (wlparam.display_coi)
    plot(glPower,freqs,'k:');
    hold on;
    plot(glPower2,freqs,'k-');
else
    plot(glPower,freqs,'k-');
end;
set(ax(3),'ylim',freqlim,'yticklabel',[],'ygrid','on');
set(ax(3),'ytick',get(ax(2),'ytick'));
if (plotstyle==1), image_reverseyaxis(ax(3)); end;

% Set current axis to scalogram before exiting
axes(ax(2));
