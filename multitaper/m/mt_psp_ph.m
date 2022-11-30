function h = mt_psp_ph(sp12,params,maxf,sp11,sp22,highlight)
%function mt_psp_ph(sp12,params,[maxf,[sp11,sp22,[highlight]]]);
%
% Phase plotting routine for use with
%   mt_sp, mt_sp2, mt_sp2_pp
%
% Plot multitaper estimates with jackknife errors if available
% Plot Neurospec confidence limits if auto-spectra provided
%
% The highlight option allows specific sub-regions of phase to be
% displayed in isolation.
%
% Input parameters
%       sp12        Cross-spectrum
%       params      Parameters structure
%       maxf        (opt) Maximum plot frequency
%       sp11        (opt) Auto-spectrum ch.1 (required for conf. limits)
%       sp22        (opt) Auto-spectrum ch.2 (      -     "     -      )
%       highlight   (opt) Vector of plot frequencies [fmin fmax]
%                           Column-vector of doublets (M x 2)
%
%function mt_psp_ph(sp12,params,[maxf,[sp11,sp22,[highlight]]]);

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

% Determine highlight vector (if empty plot whole region)
if (~exist('highlight'))
    highlight=[];
end;
if (isempty(highlight))
    highlight=[0 freqs(end)];
end;

% Phase - (NB: Neurospec conf. limits are also unbound, plot is y-limited)
for ind=1:size(highlight,1)
    % Determine highlight region
    frange2=[dsearchn(freqs',highlight(ind,1)'):dsearchn(freqs',highlight(ind,2)')];
    % Plot highlight region
    if (isfield(params,'ph'))
        if (params.jackknife)
            if (params.bootstrap)
                h(end+1)=plot(freqs(frange2),params.ph(frange2)+1.96*sqrt(params.bsphv(frange2))); hold('on');
                h(end+1)=plot(freqs(frange2),params.ph(frange2)-1.96*sqrt(params.bsphv(frange2)));
            else
                h(end+1)=plot(freqs(frange2),params.ph(frange2)+1.96*sqrt(params.jkphv(frange2))); hold('on');
                h(end+1)=plot(freqs(frange2),params.ph(frange2)-1.96*sqrt(params.jkphv(frange2)));
            end;
        end;
        plot([0 freqs(frange2(end))],[0 0],'k'); hold('on');
        ph=params.ph(frange2);
    else
        % plot([0 freqs(end)],[0 0],'k--'); hold('on');
        ph=angle(sp12(frange2)); hold('on');
        if (exist('sp22'))
            coh=abs(sp12).^2./(sp11.*sp22);
            ph95=1.96*sqrt((1./coh(frange2)-1)./(2*params.L(frange2)));
            plot(freqs(frange2),ph-ph95,'k:',freqs(frange2),ph+ph95,'k:');
        end;
    end;
    % Plot phase
    h0 = plot(freqs(frange2),ph,'k');
    % Include linear regression
    if (~((highlight(1)==0) && (highlight(2)==freqs(end))))
        % Plot regression
        regression=polyval(polyfit(freqs(frange2)',ph,1),freqs(frange2([1 end])));
        h2=plot(freqs(frange2([1 end])),regression,'k');
        set(h2,'linewidth',2);
        % Output regression angle
        disp(['Gradient (' num2str(highlight(1)) '-' num2str(highlight(2)) 'Hz): ' ...
             num2str(diff(freqs(frange2([1 end])))/diff(regression))]);
    end;
end;

% Format plot
ylim([-pi +pi]);
xlim([0 freqs(end)]);
xlabel('FREQ (Hz)');
hold('off');
title(['ph' what]);

% Format confidence limits
set(h,'color',0.6*[1 1 1]);       % Solid gray line
box('on');

h = h0;