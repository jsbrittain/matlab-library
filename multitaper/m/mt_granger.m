function I = mt_granger(sp11,sp22,sp12)

% Form spectral matrix
S = [ shiftdim(sp11,-2) shiftdim(sp12,-2); shiftdim(conj(sp12),-2) shiftdim(sp22,-2) ];

% Granger causality
I = mt_granger_spectral( S );
