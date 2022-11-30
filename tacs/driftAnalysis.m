function analysis = driftAnalysis( phase1, phase2, amp, bincount )

% Discretised phase-difference vector
phasediff = angle(exp(1i*( phase1 - phase2 )));
phasediffk = round((phasediff+pi)/2/pi*bincount);
phasediffk( phasediffk==0 ) = bincount;

% Instantaneous frequency of phase1 (accelerometer) in radians
instfreq = [ NaN; angle(exp(1i*( diff(phase1) ))) ];

% Phase bins
analysis.bincount = bincount;
analysis.phases = (1:bincount)/bincount*2*pi-pi;

% Per bin metrics
analysis.amp = nan(1,bincount);
analysis.likeli = nan(1,bincount);
analysis.instfreq = nan(1,bincount);
analysis.histcount = nan(1,bincount);
for k = (1:bincount)
    analysis.amp(k) = nanmean( amp(phasediffk==k) );
    analysis.likeli(k) = nanmean( phasediffk==k )/bincount;
    analysis.instfreq(k) = nanmedian( instfreq(phasediffk==k) );
    analysis.histcount(k) = sum( phasediffk==k );
end;

% Summary stats
analysis = driftAnalysisSummaryStats( analysis );
