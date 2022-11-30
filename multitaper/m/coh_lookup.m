function [cohmag,cohprob] = coh_lookup( cohval, K )
%function [cohmag,cohprob] = coh_lookup( cohval, K )
%
%
%
%function [cohmag,cohprob] = coh_lookup( cohval, K )

% Check for matrix input and vectorise
if (~isvector(cohval))
    % Pointwise operation so row/column-wise irrelevant
    cohmag = zeros(size(cohval)); cohprob = cohmag;
    for k = (1:size(cohval,1))
        [cohmag(k,:),cohprob(k,:)] = coh_lookup( cohval(k,:), K );
    end;
    return;
end;

lookup = load('cohlookup');
if (sum(ismember(lookup.Mlist,K))==0)
    error([' Supported taper count range: [ ' num2str(lookup.Mlist) ' ]']);
end;
if ((min(cohval)<0) || (max(cohval)>1))
    error(' Estimated Coherence range outside [0, 1]');
end;

% Convolve surface
K = find(ismember(lookup.Mlist,K));
if ( 0 )%&& exist('conv2'))
    switch ( 0 )
        case 0, kernel = 1;
            
        case 1, kernel = [0 1 0; 1 1 1; 0 1 0];
        case 2, kernel = ones( 3 );
        case 3, kernel = ones( 5, 1 );
        case 4, kernel = [ 0 0 1 0 0; 0 1 1 1 0; 1 1 1 1 1; 0 1 1 1 0; 0 0 1 0 0 ];
            
        case 11, kernel = ones( 3, 3, 3 );
    end;
    coh = convn( lookup.coh(:,:,K), kernel/sum(kernel(:)), 'same' );
else
    coh = lookup.coh(:,:,K);
end;

% Lookup (discrete precision; no interpolation)
if (size(cohval,2)>1)
    cohval = cohval.';
end;
ix0 = dsearchn( lookup.cohvals', cohval );
[cohprob, ix] = max( coh(ix0,:) , [], 2);
cohmag = lookup.cohvals(ix);

if ( 0 )
    disp( [ 'Given estimated coherence = ' num2str(cohval) ' for K = ' num2str(K) ', likely real coherence = ' num2str(cohmag) ' (p = ' num2str(cohprob) ')' ] );
end;
