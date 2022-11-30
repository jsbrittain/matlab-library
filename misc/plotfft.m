function plotfft(spect,dt,maxf);
% function plotfft(spect,[dt,[maxf]]);
%
% Plot the absolute +ve and -ve components of the spectra of a signal.
%
% Example (plot autospectra of x to 100Hz)
%       plotfft(abs(fft(x)).^2,1/rate,100);
%
% Fourier decomposition            Rearranged
% ---------------------            ----------
%
%     ||        ||                    |||    
%     |\        /|                    /|\    
%     | \      / |      ====>        / | \   
%     |  \    /  |                  /  |  \  
%     |___\--/___|                -/___|___\-
%     0         -> f            -f <-  0  -> f
%

% Assign 
if (~exist('dt'))
    dt=1;
end;
if (~exist('maxf'))
    maxf=1/dt/2;
end;

% Ensure column vector
if size(spect,1)==1
    spect=spect';
end;

% Rearrange and plot absolute spectra
n = floor((length(spect)-1)/2);
xvalues=[-n:n]/length(spect)/dt;
plot(xvalues,(abs([spect((n+2):2*n+1); spect(1:n+1)])),'k')
xlim([max([min(xvalues) -maxf]) min([max(xvalues) maxf])]);
xlabel('Hz');
%ylabel('dB');
title('Spectrum');
