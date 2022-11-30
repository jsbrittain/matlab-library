function analysis = driftAnalysisSummaryStats( analysis )

% PSI summaries for amplitude and entrainment modulation
analysis.amp_psi = mean( (analysis.amp/sum(analysis.amp)).*exp(1i*analysis.phases) );
analysis.likeli_psi = mean( (analysis.likeli/sum(analysis.likeli)).*exp(1i*analysis.phases) );
analysis.instfreq_psi = mean( (analysis.instfreq/sum(analysis.instfreq)).*exp(1i*analysis.phases) );
