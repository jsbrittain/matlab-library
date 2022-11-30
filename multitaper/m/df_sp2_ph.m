function [ph11,ph22,ph12]=df_sp2_ph(DF11,DF22,DF12)
%function [ph11,ph22,ph12]=df_sp2_ph(DF11,DF22,DF12);
%
% Dual-frequency spectrum
%
% Routine determines coherence estimates based on dual-frequency spectra
%
%function [ph11,ph22,ph12]=df_sp2_ph(DF11,DF22,DF12)

% Determine parameters
fmax=length(DF11);

% Extract stationary spectra
S11=diag(DF11);
S22=diag(DF22);

% Calculate dual-frequency coherence
allS11=S11(:,ones(1,fmax));
allS22=S22(:,ones(1,fmax));
ph11=angle(DF11);
ph22=angle(DF22);
ph12=angle(DF12);
