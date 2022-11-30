function scalogram(time,power,wlparam,opt);
%function scalogram(time,power,wlparam,[opt]);
%
% Scalogram plotting routine
%
% Required parameters
%   time    Time vector
%   power   Power of the wavelet transform
%   wlparam Wavelet parameters
%   opt     (opt) Options string
%               i        Image plot (default: contour)
%               c<style> COI plot style (0:Line 1:Fill)
%               v        Variance normalised
%
% This routine may also be used to display partial scalograms
% so long as the input parameters are accordingly scaled (time/Wsp);
% wlparam options will be scaled automatically (i.e.: wlparam.coi)
%
%function scalogram(time,power,wlparam,[opt]);

% Routine written to support frequency plotting, therefore wavelet
% transforms performed using wlTransform have frequency content added
% determined before plotting as frequencies.

% Check inputs
if (nargin<3)
    error('Not enough input parameters.');
end;
if (nargin<4)
    opt='';
end;
if (nargin>4)
    error('Too many input parameters.');
end;
if (length(time)~=size(power,2))
    error('Time vector of inconsistent length with regards wavelet transform');
end;

% Default parameters
plotstyle=0;            % 0-contour, 1-image
coistyle=0;             % 0-line, 1-fill

% Parse options string
options=deblank(opt);
opt_str='';
while (any(options))                % Parse individual options from string.
	[opt,options]=strtok(options);
	optarg=opt(2:length(opt));      % Determine option argument.
	switch (opt(1))
        case 'i'                    % Image plot
            plotstyle=1;
        case 'c'                    % COI plot style
            coistyle=str2num(optarg);
    	otherwise
            opt_str=[opt_str ' ' opt];
    end;
end;

% Determine frequency parameters if original analysis was conducted in scale
if (~isfield(wlparam,'freqs'))
    wlparam.freqs=1./(wlparam.scale*wlparam.fourier_factor);
    wlparam.coi=[Inf (1./(wlparam.fourier_factor*wlparam.coi(wlparam.coi>0))) Inf];
end;

% Determine COI plot region - locate end (Inf) points - routine may be plotting a partial scalogram
if (wlparam.display_coi)
	tstart=1;
	tstop=length(time);
	if (wlparam.coi(1)==Inf)
        tstart=2;
	end;
	if (wlparam.coi(end)==Inf)
        tstop=length(time)-1;
	end;
	tcoi=time(tstart:tstop);
	coi=wlparam.coi(wlparam.coi<Inf);    
end;

% Common parameters
mintime=min(time);
maxtime=max(time);
minfreq=min(wlparam.freqs);
maxfreq=max(wlparam.freqs);

% If required, plot on a linear scale
if (wlparam.linearscale)
    % Contour plot the wavelet transform with a Fourier frequency y-axis
    switch (plotstyle)
        case 0, contour(time,flipud(wlparam.freqs),power);
        case 1, imagesc(time,flipud(wlparam.freqs),flipud(power));
    end;
	ax(1)=gca;
	set(ax(1),'YLim',[minfreq maxfreq],'XLim',[mintime maxtime],'ygrid','on', 'xgrid','on');
	
	% Add COI
	if (wlparam.display_coi)
        % coi not provided for first and last points (default set to Inf)
        hold on;
        switch plotstyle
            case 0, plot(tcoi,coi,'k');     % Contour plot
            case 1,                         % Image plot
                switch (coistyle)
                    case 0, plot(tcoi,wlparam.freqs(end)+wlparam.freqs(1)-coi,'k'); % Line COI
                    case 1, % Shaded region COI
                        tcoi=[mintime tcoi(2:end-1) maxtime maxtime mintime mintime];
                        coi_region=wlparam.freqs(end)+wlparam.freqs(1)-[maxfreq coi(2:end-1) maxfreq minfreq minfreq maxfreq];
                        fh=fill(tcoi,coi_region,'w');
                        set(fh,'alphadatamapping','direct','facealpha',0.5);
                end;
        end;
	end;
    ylabel('Frequency (Hz)');
else
	% Contour plot the wavelet transform with a Fourier frequency y-axis
	switch (plotstyle)
        case 0, contourf(time,log2(flipud(wlparam.freqs)),power);
        case 1, imagesc(time,log2(flipud(wlparam.freqs)),flipud(power));
	end;
	ax(1)=gca;
	Yticks = 2.^(fix(log2(minfreq)):fix(log2(maxfreq)));
	set(ax(1),'YLim',log2([minfreq maxfreq]),'YTick',log2(Yticks(:)),'YTickLabel',Yticks,...
        'XLim',[mintime maxtime],'ygrid','on', 'xgrid','on');
	
	% Add COI
	if (wlparam.display_coi)
        % coi not provided for first and last points (default set to Inf)
        hold on;
        switch (plotstyle)
            case 0, plot(tcoi,log2(coi),'k');   % Contour plot
            case 1,                             % Image plot
                switch (coistyle)
                    case 0, plot(tcoi,log2(wlparam.freqs(1))+log2(wlparam.freqs(end))-log2(coi),'k'); % Line COI
                    case 1, % Shaded region COI
                        tcoi=[mintime tcoi(2:end-1) maxtime maxtime mintime mintime];
                        coi_region=log2(wlparam.freqs(1))+log2(wlparam.freqs(end))-log2([maxfreq coi(2:end-1) maxfreq minfreq minfreq maxfreq]);
                        fh=fill(tcoi,coi_region,'w');
                        set(fh,'alphadatamapping','direct','facealpha',0.5);
                end;
        end;
	end;
	
	% Add second scale axis on right-hand side corresponding to scale
	if (wlparam.display_scale)
		ax(2) = axes('position',get(ax(1),'position'));
        Yticks = 2.^(fix(log2(min(wlparam.scale))):fix(log2(max(wlparam.scale))));
		set(ax(2),'YAxisLocation','right','color','none','xgrid','off','ygrid','off','box','off','xticklabel',[], ...
		          'YLim',log2([min(wlparam.scale) max(wlparam.scale)]),'YDir','reverse','YTick',log2(Yticks(:)),'YTickLabel',Yticks,'XTick',[]);
		ylabel('Scale');
        
        % Set current handle to contour plot
        axes(ax(1));
	end;
    ylabel('Frequency (Hz)');
end;

% Adjust y axis if image plotting
if (plotstyle==1)
    image_reverseyaxis(ax(1));
end;
