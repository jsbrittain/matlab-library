function cfplot_phamp_amp( phaseangle, rmag, permrmag, amp_norm_secs, pvalue )

if (~exist('pvalue','var'))
    pvalue = [];
end;
if (isempty(pvalue))
    pvalue = 0.05;
end;
confint = icdf('norm',(1-pvalue/2),0,1);

bar(phaseangle,rmag); xlim(phaseangle([1 end]));
hold('on');
if (~isempty(amp_norm_secs))
    plot(xlim,[1 1],'k--');
end;
if (size(permrmag,1)>0)
    fill( [ phaseangle fliplr(phaseangle) ], [mean(permrmag,1)-confint*std(permrmag,[],1) fliplr(mean(permrmag,1)+confint*std(permrmag,[],1))], 'g', 'facecolor', [1 1 1]*0.8, 'facealpha', 0.8 );
end;
plot(xlim,mean(rmag)*[1 1],'r--');
title('TREMOR AMPLITUDE');
