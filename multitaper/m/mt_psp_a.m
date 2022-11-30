function ph=mt_psp_a(sp11,params,opt_str)
%function h=mt_psp_a(sp11,params,opt_str);
%
% Spectral plotting routine for use with
%   mt_sp, mt_sp2, mt_sp2_pp
%
% Plot multitaper estimates with jackknife errors if available
%
% `logplot' does not currently work with jackknife estimates
%
% Parameters
%       f<n>            Maximum frequency plot range
%       l               Linear (instead of log) plot
%       !!! b<low,high>     Blank (notch) range
%
%function h=mt_psp_a(sp11,params,opt_str);

% Default parameters
logplot = true;
maxf = [];
blankrange = [];

if (~exist('opt_str'))
    opt_str = '';
end;

% Parse options string
options=deblank(opt_str); opt_str='';
while (any(options))                % Parse individual options from string.
	[opt,options]=strtok(options);
	optarg=opt(2:length(opt));      % Determine option argument.
	switch (opt(1))
	    case 'f'                    % Specify maximum display frequency
        	n=str2num(optarg);
			if (isscalar(n))
                maxf = n;
            else
			    error(['Error in option argument -- f' optarg]);
			end;
        case 'l'                    % Linear plot axis
            logplot = false;
        case 'b'                    % Blank range
            
        otherwise                   % Options for wavelet analysis
            error(['Unknown option -- ' opt]);
    end;
end;

% Determine common parameters
if (params.jackknife)
    % Spectral confidence limits
    p=0.05;
    ttest=tq(1-p/2,params.jkcount-1);
end;
% Handles for conf. lines
h=[];

% Check if label specified
if (~isfield(params,'what'))
    what='';
else
    what=[': ' params.what];
end;

% Determine maximum plot frequency
if (isempty(maxf))
    maxf=params.rate/2;
end;
fmin=find(params.freqs>0); fmin=fmin(1);
fmax=find(min(abs(params.freqs-maxf))==abs(params.freqs-maxf));
frange=fmin:fmax;
freqs=params.freqs(fmin:fmax);

% Auto-spectra
if (params.jackknife)
    if (params.bootstrap)
        h(end+1)=plot(freqs,log10(sp11(frange))-ttest*sqrt(params.bs11v(frange))); hold('on');
        h(end+1)=plot(freqs,log10(sp11(frange))+ttest*sqrt(params.bs11v(frange)));
    else
        h(end+1)=plot(freqs,log10(sp11(frange))-ttest*sqrt(params.jk11v(frange))); hold('on');
        h(end+1)=plot(freqs,log10(sp11(frange))+ttest*sqrt(params.jk11v(frange)));
        %plot(freqs,mean(params.pv11(frange,:),2),'g-.');       % Pseudo-value estimate of log spectra
    end;
end;
% Empirical confidence limits
if (isempty(find(diff(params.L))))      % If all L the same
    xpos=freqs(1)+0.95*(freqs(end)-freqs(1))*[1 1];
    ypos=min(log10(sp11(frange)))+0.95*(max(log10(sp11(frange)))-min(log10(sp11(frange))));
    c95=2*(0.851./sqrt(params.L(1)));
    if (logplot)
        hs=plot(xpos,ypos-[0 c95],'k'); hold('on'); set(hs,'linewidth',2);
    else
        plot(freqs,10.^(log10(sp11(frange))-c95/2),'k--'); hold('on');
        plot(freqs,10.^(log10(sp11(frange))+c95/2),'k--');
    end;
end;
if (isfield(params,'bandwidth'))
    xpos=freqs(1)+0.95*(freqs(end)-freqs(1))*[1 1];
    ypos=min(log10(sp11(frange)))+0.95*(max(log10(sp11(frange)))-min(log10(sp11(frange))))-c95/2;
    plot(xpos+[-1 1]*params.bandwidth/2,ypos*[1 1],'k'); hold('on');
end;
% Plot log10 spectra
if (logplot)
    ph=plot(freqs,log10(sp11(frange)),'k');
    ylabel('x10 dB');
else
    ph=plot(freqs,sp11(frange),'k');
end;
hold('off');
xlabel('FREQ (Hz)');
axis('tight');
xlim([0 freqs(end)]);
title(['fa' what]);

% Format confidence limits
set(h,'color',[128 128 128]/255);       % Solid gray line
