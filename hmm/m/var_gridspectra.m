function var_gridspectra(hmmsp,freqs,rate,chlabels,colstr,smoother)
%function var_gridspectra(hmmsp,freqs,rate,chlabels,colstr,smoother)
%
% Grid plot spectra and coherence
%
%function var_gridspectra(hmmsp,freqs,rate,chlabels,colstr,smoother)

% Check input parameters
m=size(hmmsp,1);
if (~exist('colstr'))
    colstr = 'b';
end;
if (~exist('chlabels'))
    chlabels=[];
end;
if (isempty(chlabels))
    chlabels = num2cell(num2str(1:m));
end;
if (~exist('smoother'))
    smoother = 1;
end;

% Default parameters
logspectra = true;
maxf       = rate/2;
maxcoh     = 1;

% Plot grid-array
ylims = [Inf -Inf];
for condition = (1:size(hmmsp,4))
    col = colstr(condition);
    for ch1=(1:m)
        for ch2=(1:m)
            subplot(m,m,(ch1-1)*m+ch2); hold('on');
            if (ch1==ch2)       % Diagonal (spectra)
                if (logspectra)
                    sp = smooth(10*log10(squeeze(hmmsp(ch1,ch1,:,condition))),smoother);
                else
                    sp = smooth(squeeze(hmmsp(ch1,ch1,:,condition)),smoother);
                end;
                plot(freqs,sp,col);
                ylims = [min([sp; ylims(1)]) max([sp; ylims(2)])];
                xlim([0 maxf]); title([num2str(ch1) ': ' chlabels{ch1}]);
            elseif (ch1<ch2)    % Upper quadrant (Coherence)
                plot(freqs,smooth(squeeze(abs(hmmsp(ch1,ch2,:,condition)).^2./(hmmsp(ch1,ch1,:,condition).*hmmsp(ch2,ch2,:,condition))),smoother),col);
                box('on'); xlim([0 maxf]); ylim([0 maxcoh]); title(''); xlabel('');
            else                % Lower quandrant (cross-spectra  / Phase)
                switch (1)
                    case 0,
                        plot(freqs,10*log10(smooth(squeeze(abs(hmmsp(ch1,ch2,:,condition))),smoother)),col);
                    case 1,
                        plot(freqs,smooth((angle(squeeze(hmmsp(ch1,ch2,:,condition)))),smoother),col);
                        ylim([-pi pi]);
                end;
                title(''); xlabel(''); xlim([0 maxf]);
            end;
        end;
    end;
end;

for ch1=(1:m)
    subplot(m,m,(ch1-1)*m+ch1);
    ylim(ylims);
end;
