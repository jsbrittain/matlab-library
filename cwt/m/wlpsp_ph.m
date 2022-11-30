function wlpsp_ph(time,Ph,wlparam,opt_str);
%function wlpsp_ph(time,Ph,wlparam,[opt_str]);
%
% Plot wavelet phase as an image
%
% Input parameters
%   time        Time vector
%   Ph          Wavelet phase
%   wlparam     Wavelet transform parameters
%   opt_str     (opt) Options string
%                   c - Display colour bar (default: no)
%
%function wlpsp_ph(time,Ph,wlparam,[opt_str]);

% Default parameters
cb_display=0;

% Check input parameters
if (exist('cl'))
    cl_limits=1;
    if (length(cl)>1)
        cl_window=1;
	end;
end;
if (nargin<4)
    opt_str='';
end;

% Parse options string
options=deblank(opt_str);
while (any(options))                % Parse individual options from string.
	[opt,options]=strtok(options);
	optarg=opt(2:length(opt));      % Determine option argument.
	switch (opt(1))
        case 'c'
            cb_display=1;
    	otherwise
            warning([' (wlpsp_ph.m) Unknown option: ' opt]);
    end;
end;

% Get colomap properties
cmap=colormap;
k=size(cmap,1);

% Plot phase
if (wlparam.linearscale)
    imagesc(time([1 end]),wlparam.freqs([1 end]),flipud(Ph));
else
    imagesc(time([1 end]),log2(wlparam.freqs([1 end])),flipud(Ph));
    % Determine plot parameters
    mintime=min(time);
	maxtime=max(time);
    minfreq=min(wlparam.freqs);
    maxfreq=max(wlparam.freqs);
    % Define logarithmic plot axis
    ax(1)=gca;
	Yticks = 2.^(fix(log2(minfreq)):fix(log2(maxfreq)));
	set(ax(1),'YLim',log2([minfreq maxfreq]),'YTick',log2(Yticks(:)),'YTickLabel',Yticks,'XLim',[mintime maxtime]);
end;

% Label plot
image_reverseyaxis;
title('Phase');

% Add colour bar
if (cb_display)
	ch=colorbar;
    set(ch,'ytick',[-pi -pi/2 0 pi/2 pi]);
	set(ch,'yticklabel',{'-pi','-pi/2',' 0',' pi/2',' pi'});
    % Rescale axes to reduce colorbar influence
    ah=gca;
    ahpos=get(ah,'position');
    chpos=get(ch,'position');
    width=1/4;                      % Rescale to a portion of original colorbar width
    ahpos(3)=ahpos(3)+(1-width)*chpos(3)+(chpos(1)-ahpos(1)-ahpos(3))/2;
    chpos(1)=chpos(1)+(1-width)*chpos(3);
    chpos(3)=width*chpos(3);
    set(ah,'position',ahpos);
    set(ch,'position',chpos);
end;

% Add COI
if (wlparam.display_coi)
    hold on;
    tcoi=time(wlparam.coi<Inf);
    coi=wlparam.coi(wlparam.coi<Inf);
    if (wlparam.linearscale)
        plot(tcoi,wlparam.freqs(end)+wlparam.freqs(1)-coi,'k');
    else
        plot(tcoi,log2(wlparam.freqs(1))+log2(wlparam.freqs(end))-log2(coi),'k');
    end;
end;
