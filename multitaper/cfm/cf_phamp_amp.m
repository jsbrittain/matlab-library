function cf_phamp_amp( phaseangle, rmag, permrmag, amp_norm_secs )

bar(phaseangle,rmag); xlim(phaseangle([1 end]));
hold('on');
if (~isempty(amp_norm_secs))
    plot(xlim,[1 1],'k--');
end;
if (size(permrmag,1)>0)
    fill( [ phaseangle fliplr(phaseangle) ], [mean(params.perm.rmag,1)-2*std(params.perm.rmag,[],1) fliplr(mean(params.perm.rmag,1)+2*std(params.perm.rmag,[],1))], 'g', 'facecolor', [1 1 1]*0.8, 'facealpha', 0.8 );
end;
plot(xlim,mean(params.rmag)*[1 1],'r--');
title('TREMOR AMPLITUDE');
