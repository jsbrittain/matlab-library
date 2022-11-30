function wlpsp_coh(time,WCo,wlparam,arg1,arg2)
%function wlpsp_coh(time,WCo,wlparam,[cl],[opt_str])
%
% Plot wavelet coherence as an image
%
% Input parameters
%   time        Time vector
%   WCo         Wavelet coherence
%   wlparam     Wavelet transform parameters
%   cl          Confidence limits structure (2-cols)
%                   1: Frequency (Hz)
%                   2: 95% confidence level (based on a null hypothesis)
%   opt_str     (opt) Options string
%                   c               Display colour bar (default: no)
%                   s<cmin,cmax>    Scaling for colorbar
%
%function wlpsp_coh(time,WCo,wlparam,[cl],[opt_str])

% Default parameters
cl_limits=0;        % Confidence limits provided
cb_display=0;
clim=[];
opt_str='';

% Check input parameters
if (nargin==5)
    cl=arg1;
    opt_str=arg2;
end;
if (nargin==4)
    if (ischar(arg1))
        opt_str=arg1;
    else
        cl=arg1;
    end;
end;
if (exist('cl','var'))
    if (~isempty(cl))
        cl_limits=1;
    else
        if (isfield(wlparam,'L'))
            cl = (1-0.05^(1/(wlparam.L)));
        else
            cl = (1-0.05^(1/(wlparam.trigcount)));
        end;
        cl_limits=1;
    end;
end;

% Check input data
if (~isfield(wlparam,'freqs'))
    wlparam.freqs=1./(wlparam.fourier_factor*wlparam.scale);
    wlparam.coi=[Inf (1./(wlparam.fourier_factor*wlparam.coi(wlparam.coi>0))) Inf];
end;

% Parse options string
options=deblank(opt_str);
while (any(options))                % Parse individual options from string.
	[opt,options]=strtok(options);
	optarg=opt(2:length(opt));      % Determine option argument.
	switch (opt(1))
        case 'c'
            cb_display=1;
        case 's'                    % Scaling for colorbar [cmin,cmax]
            [clim1,clim2]=strtok(optarg,',');
            clim=[str2num(clim1) str2num(clim2(2:end))];
            clear('clim1','clim2');
        otherwise
            warning([' (wlpsp_coh.m) Unknown option: ' opt]);
    end;
end;

% Determine colormap parameters
cmap=colormap;
k=size(cmap,1);     % The number of color indices

% Specify minimum and maximum values for scaling
if (isempty(clim))    % Used for scaling normalisation
    clim=[0 max(WCo(:))];
end;

% Display coherence as an image
if (wlparam.linearscale)
    imagesc(time([1 end]),wlparam.freqs([1 end]),flipud(WCo),clim);
else
    imagesc(time([1 end]),log2(wlparam.freqs([1 end])),flipud(WCo),clim);
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

% Determine confidence limits
if (length(cl)==1)
    % Scalar confidence limits provided
    clR95=cl; clear('cl');
    cl.R95=clR95; clear('clR95');
    conf95=(WCo>=cl.R95);
elseif (~isempty(cl))
    % Convert cols to structure array cl.(freqs,R95)
    clfreqs=cl; clear('cl');
    cl.freqs=clfreqs; clear('clfreqs');
    cl.R95=cl.freqs(:,2);
    cl.freqs=cl.freqs(:,1);
    % Determine confidence limits
    conf95=false(length(wlparam.freqs),size(WCo,2));
    for ind=1:length(cl.freqs)
        conf95(ind,:)=(WCo(ind,:)>=cl.R95(ind));
    end;
    % % Set conf limit to logical(0) inside coi
    % for ind=1:length(wlparam.freqs)
    %     conf95(ind,(wlparam.freqs(ind)<wlparam.coi))=0;
    % end;
end;

% Add colorbar
if (cb_display)
    colorbar;
%     ah=gca;
%     ch=colorbar;
%     if (cl_limits)  % Draw horizontal line marking static conf. limit
%         if (isempty(find(cl.R95(1)~=cl.R95, 1)))
%             subplot(ch);
%             hold on;
%             xlims=xlim;
%             plot(xlims,[1 1]*cl.R95(1),'k');
%             subplot(ah);
%         end;
%     end;
%     % Rescale axes to reduce colorbar influence
%     ahpos=get(ah,'position');
%     chpos=get(ch,'position');
%     width=1/4;                      % Rescale to a portion of original colorbar width
%     ahpos(3)=ahpos(3)+(1-width)*chpos(3)+(chpos(1)-ahpos(1)-ahpos(3))/2;
%     chpos(1)=chpos(1)+(1-width)*chpos(3);
%     chpos(3)=width*chpos(3);
%     set(ah,'position',ahpos);
%     set(ch,'position',chpos);
end;

% Overlay confidence limits on coherence figure
if (0)%cl_limits)
    hold on;
    % Plot confidence limits (outside the coi only)
    if (length(time)==size(conf95,2))
        if (wlparam.linearscale)
            contour(time,wlparam.freqs,flipud(conf95),1,'k');
        else
            contour(time,log2(wlparam.freqs),flipud(conf95),1,'k');
        end;
    end;
end;

% Format plots
image_reverseyaxis;
title('Wavelet Coherence');
ylabel('Frequency (Hz)');

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
    if (isfield(wlparam,'plotlimit'))
        if (isfield(wlparam.plotlimit,'minfreq'))
            plot(time([1 end]),log2(wlparam.freqs(1))+log2(wlparam.freqs(end))-log2([1 1]*wlparam.plotlimit.minfreq),'k');
            plot(time([1 end]),log2(wlparam.freqs(1))+log2(wlparam.freqs(end))-log2([1 1]*wlparam.plotlimit.maxfreq),'k');
        end;
    end;
end;
