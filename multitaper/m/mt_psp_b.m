function mt_psp_b(sp22,params,maxf,logplot);
%function mt_psp_b(sp22,params,[maxf,[logplot]]);
%
% Spectral plotting routine for use with
%   mt_sp, mt_sp2, mt_sp2_pp
%
% Plot multitaper estimates with jackknife errors if available
%
% `logplot' does not currently work with jackknife estimates
%
%function mt_psp_b(sp22,params,[maxf,[logplot]]);

if (~exist('logplot'))
    logplot=true;
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
if (~exist('maxf'))
    maxf=params.rate/2;
end;
fmin=find(params.freqs>0); fmin=fmin(1);
fmax=find(min(abs(params.freqs-maxf))==abs(params.freqs-maxf));
frange=fmin:fmax;
freqs=params.freqs(fmin:fmax);

% Auto-spectra
if (params.jackknife)
    if (params.bootstrap)
        h(end+1)=plot(freqs,log10(sp22(frange))-ttest*sqrt(params.bs11v(frange))); hold('on');
        h(end+1)=plot(freqs,log10(sp22(frange))+ttest*sqrt(params.bs11v(frange)));
    else
        h(end+1)=plot(freqs,log10(sp22(frange))-ttest*sqrt(params.jk11v(frange))); hold('on');
        h(end+1)=plot(freqs,log10(sp22(frange))+ttest*sqrt(params.jk11v(frange)));
        %plot(freqs,mean(params.pv11(frange,:),2),'g-.');       % Pseudo-value estimate of log spectra
    end;
end;
% Empirical confidence limits
if (isempty(find(diff(params.L))))      % If all L the same
    xpos=freqs(1)+0.95*(freqs(end)-freqs(1))*[1 1];
    ypos=min(log10(sp22(frange)))+0.95*(max(log10(sp22(frange)))-min(log10(sp22(frange))));
    c95=2*(0.851./sqrt(params.L(1)));
    if (logplot)
        hs=plot(xpos,ypos-[0 c95],'k'); hold('on'); set(hs,'linewidth',2);
    else
        plot(freqs,10.^(log10(sp22(frange))-c95/2),'k--'); hold('on');
        plot(freqs,10.^(log10(sp22(frange))+c95/2),'k--');
    end;
end;
if (isfield(params,'bandwidth'))
    xpos=freqs(1)+0.95*(freqs(end)-freqs(1))*[1 1];
    ypos=min(log10(sp22(frange)))+0.95*(max(log10(sp22(frange)))-min(log10(sp22(frange))))-c95/2;
    plot(xpos+[-1 1]*params.bandwidth/2,ypos*[1 1],'k'); hold('on');
end;
% Plot log10 spectra
if (logplot)
    plot(freqs,log10(sp22(frange)),'k');
    ylabel('x10 dB');
else
    plot(freqs,sp22(frange),'k');
end;
hold('off');
xlabel('FREQ (Hz)');
axis('tight');
xlim([0 freqs(end)]);
title(['fa' what]);

% Format confidence limits
set(h,'color',[128 128 128]/255);       % Solid gray line
