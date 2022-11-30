function wl_psp_grid( wltime, Wsp, wlparam, baseline, chlabels )
%function wl_psp_grid( wltime, Wsp, wlparam, baseline, chlabels )
%
% Wavelet (spectrogram) grid plot
% Percetnage change auto-spectra
%

% Check time input and reform as appropriate
M = length(Wsp);
if (length(wltime)~=size(Wsp{1,1},2))
    wltime = wl_time(wltime, Wsp{1,1});
end;

% Check channel labels
if (~exist('chlabels'))
    for ind=(1:M)
        chlabels{ind}=num2str(ind);
    end;
end;

% Plot grid-array
figure;
for ch1=(1:M)
    for ch2=(1:M)
        subplot(M,M,(ch1-1)*M+ch2);
        if (ch1==ch2)       % Diagonal (spectra)
            if (any(isnan(Wsp{ch1,ch2}(:,:,1))))
                continue;
            end
            wlpsp_a(wltime, Wsp{ch1,ch2}(:,:,1), wlparam{ch1,ch2},[],'l');
            title(chlabels{ch1});
        elseif (ch1<ch2)    % Upper quadrant (Coherence)
            if (any(isnan(Wsp{ch1,ch2}(:,:,4))))
                continue;
            end
            wlpsp_coh(wltime, Wsp{ch1,ch2}(:,:,4), wlparam{ch1,ch2},[],'');
            %delete(gca);
            xlabel(''); ylabel(''); title('');
        else                % Lower quandrant (Phase)
            delete(gca);
        end;
    end;
end;
