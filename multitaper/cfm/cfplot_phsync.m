function cfplot_phsync( params )

figure;

% Boxplot / histogram of phase synchronisation index (PSI)
subplot(3,3,3);
    switch (1)
        case 0,     % Boxplot
            boxplot(params.ps);
        case 1,     % Histogram
            [y,x] = hist(params.ps,(0:0.05:1));
            barh(x,y); ylim([0 1]); hold('on');
    end;

% Colour-code plots by PSI
splitps = true;
if (splitps)
    switch (1)
        case 0,     % Constant threshold
            th_pssplit = 0.5;
        case 1,     % Median
            th_pssplit = median(params.ps);
    end;
    % Format second half of histogram plot
    y(x<th_pssplit) = nan; barh(x,y,'r');
end;


% Mean phase synchronisation index (PSI) within blocks
subplot(3,3,[1 2]);
    plot( params.time, params.ps, 'b' );
    if (splitps)
        hold('on'); pdps = params.ps; pdps(params.ps<th_pssplit) = nan;
        plot( params.time, pdps, 'r' );
    end;
    xlim(params.time([1 end])); ylim([0 1]);
    title('PHASE SYNCHRONISATION INDEX');

% Mean phase difference within blocks
subplot(3,3,[4 5]);
    plot( params.time, params.pd, 'b' );
    if (splitps)
        hold('on'); pdps = params.pd; pdps(params.ps<th_pssplit) = nan;
        plot( params.time, pdps, 'r' );
    end;
    xlim(params.time([1 end])); ylim([-pi pi]);
    title('MEAN PHASE DIFFERENCE');
    
% Mean tremor amplitude within blocks
subplot(3,3,[7 8]);
    plot( params.time, params.xmag, 'b' );
    if (splitps)
        hold('on'); pdps = params.xmag; pdps(params.ps<th_pssplit) = nan;
        plot( params.time, pdps, 'r' );
    end;
    xlim(params.time([1 end]));
    title('TREMOR MAGNITUDE');

% ( amplitude - PSI ) correlation
subplot(3,3,9);
    linefit(params.xmag,params.ps);
    [R,P]=corrcoef(params.xmag',params.ps');
    ylim([0 1]); ylabel('PSI'); xlabel('|Ax|');
    title(['R=' num2str(R(1,2)) ', P=' num2str(P(1,2))]);
