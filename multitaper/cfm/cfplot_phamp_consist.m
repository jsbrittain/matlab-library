function cfplot_phamp_consist( preprocdata, stimfreq, strTitle )
%
%
% Most preprocessing and amplitude analysis done on 'x', whereas 'y'
% provides the phase comparitor.
%
%

% Expand proprocessed data to (local) workspace
struct2vars( preprocdata );

% Reset severe deflections to mean and smooth
%ifxh(ifxh>(2*passbandx(2))) = mean(ifxh);
%ifxh = smooth(ifxh,smoothing);

%% Analyse TACS data, then randomly permute phase and continue re-analysing

figure; K = 6;
time = (0:length(ifxh)-1)/60/rate;
for n = (1:2)
    if (n==1)
        subplot(1,K,(1:K-1));
    else
        subplot(1,K,K);
    end;
    h = plotyy(time,ifxh,time,magxh);
    set(h,'xlim',time([1 end]));
    set(h(1),'ylim',[0 15]); set(h(2),'ylim',[0 1]);
    set(h(1),'ytick',[0:15]); set(h(2),'ytick',(0:0.2:1));
    axes(h(1)); hold('on');
    plot(time,smooth(ifxh,2*rate),'c');
    if (~isempty(stimfreq))
        plot(time([1 end]),stimfreq*[1 1],'r','linewidth',2);
        if (stimfreq>8)
            plot(time([1 end]),stimfreq/2*[1 1],'r--','linewidth',2);
        end;
    end;
    axes(h(2));
    if (n==1)
        title( strTitle );
        ylabel(h(1),'INSTANTANEOUS FREQUENCY');
    else
        ylabel(h(2),'TREMOR AMPLITUDE');
    end;
    xlabel('TIME (MINS)');
end;
