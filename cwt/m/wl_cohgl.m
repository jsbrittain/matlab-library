function glPower=wl_cohgl(Wsp,wlparam);
%function glPower=wl_cohgl(Wsp,wlparam);
%
% Function to generate global coherence values based on auto and cross-spectral
% estimates as provided by the Wavelet spectral coefficient matrix Wsp.
% The resultant global power neglects the region of COI if wlparam.display_coi is set.
%
%function glPower=wl_cohgl(Wsp,wlparam);

% Calculate global coherence from global auto/cross-spectra
N=size(Wsp,2);
a11=Wsp(:,:,1); a22=Wsp(:,:,2); a12=Wsp(:,:,3);
if (wlparam.display_coi)    % Remove COI region if required
	for ind=1:N
        a11(wlparam.freqs<wlparam.coi(ind),ind)=0;
        a22(wlparam.freqs<wlparam.coi(ind),ind)=0;
        a12(wlparam.freqs<wlparam.coi(ind),ind)=0;
	end;
end;
a11=mean(a11,2); a22=mean(a22,2); a12=mean(a12,2);

% Generate global coherence - Allow division by zero resulting in NaN: Such numbers are not displayed on plots
warning off
glPower=(abs(a12).^2)./(a11.*a22);
warning on
