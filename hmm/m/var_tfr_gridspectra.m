function var_gridspectra(hmmsp,freqs,rate,chlabels)
%function var_gridspectra(hmmsp,freqs,rate,chlabels)
%
% Grid plot spectra and coherence
%
%function var_gridspectra(hmmsp,freqs,rate,chlabels)

% Check input parameters
m=size(hmmsp,1);
if (~exist('chlabels'))
    chlabels=num2cell(num2str(1:m));
end;

% Default parameters
smoother   = 1;
logspectra = true;
maxf       = rate/2;
maxcoh     = 1;

% Plot grid-array
figure; clims = [Inf -Inf];
time=(1:size(hmmsp,4))/rate;
for ch1=(1:m)
    for ch2=(1:m)
        subplot(m,m,(ch1-1)*m+ch2); hold('on');
        if (ch1==ch2)       % Diagonal (spectra)
            if (logspectra)
                sp = 10*log10(squeeze(hmmsp(ch1,ch1,:,:)));
            else
                sp = squeeze(hmmsp(ch1,ch1,:,:));
            end;
            imagesc(time,freqs,sp);
            xlim(time([1 end])); ylim([0 freqs(end)]);
            title([num2str(ch1) ': ' chlabels{ch1}]);
            clims = [min([sp(:); clims(1)]) max([sp(:); clims(2)])];
        elseif (ch1<ch2)    % Upper quadrant (Coherence)
            imagesc(time,freqs,squeeze(abs(hmmsp(ch1,ch2,:,:)).^2./(hmmsp(ch1,ch1,:,:).*hmmsp(ch2,ch2,:,:))));
            box('on'); ylim([0 maxcoh]); title(''); xlabel('');
            xlim(time([1 end])); ylim([0 freqs(end)]);
        else                % Lower quandrant (cross-spectra  / Phase)
            switch (1)
                case 0,
                    imagesc(time,freqs,10*log10(squeeze(abs(hmmsp(ch1,ch2,:,:)))));
                case 1,
                    imagesc(time,freqs,(angle(squeeze(hmmsp(ch1,ch2,:,:)))));
            end;
            xlim(time([1 end])); ylim([0 freqs(end)]); title(''); xlabel('');
        end;
        %rescale_colorbar;
    end;
end;
for ch1=(1:m)
    subplot(m,m,(ch1-1)*m+ch1);
    set(gca,'clim',clims);
end;
