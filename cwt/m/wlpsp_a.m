function wlpsp_a(time,S,wlparam,arg1,arg2)
%function wlpsp_a(time,S,wlparam,[cl],[opt_str]);
%
% Display log10 autospectra with optional confidence limits
%
% Input parameters
%   time        Time vector
%   S           Autospectra
%   wlparam     Wavelet transform parameters
%   cl          (opt) Confidence limits (log10(mean(S))+/-cl) (PBMB 6.2)
%                   Scalar or 2-col vector [freq,cl]
%   opt_str     (opt) Options string
%                   g               Local significane (default: local)
%                   b<mint,maxt>    Baseline spectra between time-limits
%                   c               Include a colour bar (default: no)
%                   s<cmin,cmax>    Scaling for colorbar
%                   l               Non-logscale
%
%function wlpsp_a(time,S,wlparam,[cl],[opt_str]);

% Define default parameter options
cl_limits=false;
cb_display=false;
global_sig=false;
logscale=true;
clim=[];
opt_str='';
baseline=[];

% Determine input parameters
if (nargin==5)
    cl=arg1;
    opt_str=arg2;
end;
if (nargin==4)
    if (isstr(arg1))
        opt_str=arg1;
    else
        cl=arg1;
    end;
end;
if (exist('cl'))
    if (~isempty(cl))
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
        case 'g'                    % Fixed integration window
            global_sig=1;
        case 's'                    % Scaling for colorbar [cmin,cmax]
            [clim1,clim2]=strtok(optarg,',');
            clim=[str2num(clim1) str2num(clim2(2:end))];
            clear('clim1','clim2');
        case 'b'                    % Baseline
            [blim1,blim2]=strtok(optarg,',');
            baseline=[str2num(blim1) str2num(blim2(2:end))];
            clear('blim1','blim2');
        case 'l'                    % Non-logscale
            logscale=false;
    	otherwise
            error([' (wlpsp_a.m) Unknown option: ' opt]);
    end;
end;

% Determine colour map parameters
cmap=colormap;
k=size(cmap,1);

% Baseline spectra
if (~isempty(baseline))
    
    % Match time vector to size of spectrogram matrix
    time = time(1) + ((0:(size(S,2)-1))/(size(S,2)-1))*(time(end)-time(1));
    
    base=dsearchn(time',baseline');
    base=(base(1):base(end));
    if (logscale)
        S=S./repmat(mean(S(:,base),2),1,size(S,2));
    else
        S=S-repmat(mean(S(:,base),2),1,size(S,2));
    end;
end;

% Autospectra ch.1 confidence limit (Place contour around significant variations from the mean)
if (cl_limits)
    
    if (~logscale)
        error(' Conf lims defined for log scale only at present.');
    end;
        
    % Determine coi region
%     outsidecoi=logical(zeros(size(S)));
%     for ind=1:length(wlparam.freqs)
%         outsidecoi(ind,:)=(wlparam.freqs(ind)>=wlparam.coi);
%     end;
    outsidecoi = true(size(S));      % Use whole transform for stats (as with Fourier)
    
    % Determine time averages for various frequencies
    meanS=zeros(size(S,1),1);
    for ind=1:length(wlparam.freqs)
        meanS(ind)=mean(S(ind,outsidecoi(ind,:)));
    end;
    
    % Change to global significance
    if (global_sig)
        meanS(:)=mean(meanS);
    end;
    
    % Convert spectra (and mean/asymptotic value) to log10
    meanS=log10(meanS);
    warning off
	S=log10(S);
	warning on
    
    % Construct significance contour
    sigS= false(size(S));
    for ind=1:length(wlparam.freqs)
        % 95% confidence interval
        if (length(cl)>1)
            WSsig=cl(wlparam.freqs(ind)==cl(:,1),2);
        else
            WSsig=cl;
        end;
        
        % Determine areas outside of confidence interval
        sigS(ind,outsidecoi(ind,:))=abs(S(ind,outsidecoi(ind,:))-meanS(ind))>WSsig;
    end;
else
    % Convert spectra to log10
    if (logscale)
        warning off
        S=log10(S);
        warning on
    end;
end;

% Specify minimum and maximum values for scaling
if (isempty(clim))    % Used for scaling normalisation
    clim=[min(min(S)) max(max(S))];
end;

% Display as an image (normalised values stretched over colormap indices)
if (wlparam.linearscale)
    imagesc(time([1 end]),wlparam.freqs([1 end]),flipud(S),clim);
    % Place a contour of significant variations from the mean over the image
	if (cl_limits)
        hold on
        warning off MATLAB:conversionToLogical
        ch=contour(time,wlparam.freqs,flipud(sigS),1,'k');
        warning on MATLAB:conversionToLogical
	end;
else
    imagesc(time([1 end]),log2(wlparam.freqs([1 end])),flipud(S),clim);
    % Determine plot parameters
    mintime=min(time);
	maxtime=max(time);
    minfreq=min(wlparam.freqs);
    maxfreq=max(wlparam.freqs);
    % Define logarithmic plot axis
    ax(1)=gca;
	Yticks = 2.^(fix(log2(minfreq)):fix(log2(maxfreq)));
	set(ax(1),'YLim',log2([minfreq maxfreq]),'YTick',log2(Yticks(:)),'YTickLabel',Yticks,'XLim',[mintime maxtime]);
    % Place a contour of significant variations from the mean over the image
	if (cl_limits)
        hold on
        warning off MATLAB:conversionToLogical
        ch=contour(time,log2(wlparam.freqs),flipud(sigS),1,'k');
        warning on MATLAB:conversionToLogical
	end;
end;
image_reverseyaxis;

% Label plot and axes
title('Autospectra');
ylabel('Frequency (Hz)');

% Add colour bar
if (cb_display)
    colorbar;
%     % Create colorbar
%     ah=gca;
% 	ch=colorbar;
%     subplot(ch);
%     title('x10 dB');
%     subplot(ah);
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

% Report percentage area covered by signifcant fluctuations from the mean
%disp(['Area covered by 95% significant fluctuations from the mean: ' num2str(100*length(sigS(sigS==1))/length(outsidecoi(outsidecoi==1))) '%']);
