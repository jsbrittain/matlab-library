function analysis = phaseprecession( acc, tacs, bincount )

hacc = hilbert(acc);
htacs = hilbert(tacs);
phasediff = angle(exp(1i*(angle(hacc)-angle(htacs))));
magacc = abs(hacc);

phasediff = round( bincount*(phasediff+pi)/2/pi );
phasediff(phasediff==0) = bincount;

analysis = [];
analysis.likeli = zeros(1,bincount);
analysis.driftamp = zeros(1,bincount);
analysis.driftmed = zeros(1,bincount);

for k = (1:bincount)
    analysis.likeli(k)   = sum( (phasediff==k) )/length(phasediff);
    analysis.driftamp(k) = mean( magacc(phasediff==k) );
    analysis.driftmed(k) = median( magacc(phasediff==k) );
end;

analysis.phases = (2*pi*(1:bincount)/bincount-pi) - pi/bincount;     % Extra -pi/bincount centralises phases

if (nargout==0)
    figure;
    subplot(2,3,1);
        bar( analysis.phases, analysis.likeli );
        title('Likelihood'); xlabel('{\it\phi}_{diff}');
    subplot(2,3,2);
        bar( analysis.phases, analysis.driftamp );
        title('Mean amplitude'); xlabel('{\it\phi}_{diff}');
    subplot(2,3,3);
        bar( analysis.phases, analysis.driftmed );
        title('Median amplitude'); xlabel('{\it\phi}_{diff}');
end;
