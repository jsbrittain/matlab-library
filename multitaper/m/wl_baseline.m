function S = wl_baseline( time, S, timerange )
%function S = wl_baseline( time, S, timerange )
%
% Percentage change from baseline
%
%
%function S = wl_baseline( time, S, timerange )

trange = dsearchn( time.', timerange.' );

for ch = (1:size(S,3))
    basespectrum = repmat( mean( S(:,trange(1):trange(2),ch), 2 ), 1, size(S,2) );
    S(:,:,ch) = (S(:,:,ch) - basespectrum)./basespectrum;
end;