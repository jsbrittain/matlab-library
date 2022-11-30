function jk_comp_pcoh(coh1,coh2,dz,dzv,params);
%
% Function to plot coherence comparison from jk_comp_coh.m
%

% Coherence pair
figure;
subplot(3,1,1);
plot(params.freqs,coh1,'k',params.freqs,coh2,'b');
legend('Pair 1','Pair 2');
title('Coherence');

% Normalised difference of coherence
subplot(3,1,2);
plot(params.freqs,dz,'k');
hold on;
plot(params.freqs(abs(dz)>1.96),dz(abs(dz)>1.96),'k.');
plot(params.freqs([1 end]),[0 0],'k');
plot(params.freqs([1 end]),1.96*[1 1],'k:',params.freqs([1 end]),-1.96*[1 1],'k:');
legend('\Deltaz','Reject H0');
ylabel('\Deltaz');

% Variance
subplot(3,1,3);
plot(params.freqs,dzv,'k',params.freqs([1 end]),[1 1],'k');
ylabel('\sigma_J^2');
xlabel('Freq (Hz)');
