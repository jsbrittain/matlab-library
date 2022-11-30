function cfplot_phamp_amp_instfreq( phaseangle, rif, permrif, pvalue )

if (~exist('pvalue','var'))
    pvalue = [];
end;
if (isempty(pvalue))
    pvalue = 0.05;
end;
confint = icdf('norm',(1-pvalue/2),0,1);

bar(phaseangle,rif); hold('on');
fill( [ phaseangle fliplr(phaseangle) ], [mean(permrif,1)-confint*std(permrif,[],1) fliplr(mean(permrif,1)+confint*std(permrif,[],1))], 'g', 'facecolor', [1 1 1]*0.8, 'facealpha', 0.8 );
xlim([-pi pi]); %ylim([min(rif) max(rif)]);
plot(xlim,mean(rif)*[1 1],'r--');
title('INSTANTANEOUS FREQUENCY');
