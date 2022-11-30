function [coh11,coh22,coh12]=df_sp2_coh(DF11,DF22,DF12);
%function [coh11,coh22,coh12]=df_sp2_coh(DF11,DF22,DF12);
%
% Dual-frequency spectrum
%
% Routine determines coherence estimates based on dual-frequency spectra
%
%function [coh11,coh22,coh12]=df_sp2_coh(DF11,DF22,DF12);

% Determine parameters
fmax=length(DF11);

% Extract stationary spectra
S11=diag(DF11);
S22=diag(DF22);

% Calculate dual-frequency coherence
allS11=S11(:,ones(1,fmax));
allS22=S22(:,ones(1,fmax));
coh11=abs(DF11).^2./(allS11.*allS11');
coh22=abs(DF22).^2./(allS22.*allS22');
coh12=abs(DF12).^2./(allS22.*allS11');
