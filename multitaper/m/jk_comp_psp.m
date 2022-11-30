function jk_comp_psp(spy11,spy22,dy,dyv,params);
%
% Function to plot spectral comparison from jk_comp_sp.m
%

% Plot results
figure;
subplot(3,1,1);
plot(params.freqs,spy11,'k',params.freqs,spy22,'b');
ylabel('dB');
legend('Ch.1','Ch.2');
title('Spectra');

% Normalised difference of spectra
subplot(3,1,2);
plot(params.freqs,dy,'k');
hold on;
plot(params.freqs(abs(dy)>1.96),dy(abs(dy)>1.96),'k.');
plot(params.freqs([1 end]),[0 0],'k');
plot(params.freqs([1 end]),1.96*[1 1],'k:',params.freqs([1 end]),-1.96*[1 1],'k:');
legend('\Deltay','Reject H0');
ylabel('\Deltay');

% Variance
subplot(3,1,3);
plot(params.freqs,dyv,'k',params.freqs([1 end]),[1 1],'k');
ylabel('\sigma_J^2');
xlabel('Freq (Hz)');
