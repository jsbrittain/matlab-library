function cfplot_phamp_amp_incidence( phaseangle, rcount, permrcount, pvalue )

if (~exist('pvalue','var'))
    pvalue = [];
end;
if (isempty(pvalue))
    pvalue = 0.05;
end;
confint = icdf('norm',(1-pvalue/2),0,1);

bar(phaseangle,rcount);
hold('on'); plot(phaseangle,mean(rcount)); xlim([-pi pi]);
if (size(permrcount,1)>0)
    fill( [ phaseangle fliplr(phaseangle) ], [ mean(permrcount,1)-confint*std(permrcount,[],1) fliplr(mean(permrcount,1)+confint*std(permrcount,[],1))], 'g', 'facecolor', [1 1 1]*0.8, 'facealpha', 0.8 );
end;
title('INCIDENCE');