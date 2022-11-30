
% NB: The number of oscillators (say G) produces a triangular number of
% phase-difference pairs INdim. These are mapped onto a M = Indim + 1
% Cartesian plane, where the extra-dimension comes from the resitrcition of
% amplitude to r = 1.

Nnodes = 3;
INdim = Nnodes*(Nnodes-1)/2;        % Number of pairs of oscillators (triangular number)

N = 1e3;
coupling = 0.1;                     % This is not (necessarily) a real coupling strength!
ph = 2*pi*rand(N,INdim) + repmat(2*pi*rand(1,INdim),N,1);
subdim = (1:size(ph,2));
ph(:,subdim) = ( 1 - coupling )*2*pi*rand( N, length(subdim) ) + repmat(2*pi*rand(1,length(subdim)),N,1);

% Calculate { n-dimensional PSI, n-dim average, n-dim elements }
[psi,psid,x] = psin( ph );

% Plot first 2-dimensions
figure;
subplot(121);
    plotn( x, '.' ); plotn( [zeros(size(psid)); psid], 'r.-', 'markersize', 10 );
    title(sprintf('%gD plot of %g-channel-difference phase data, PSI = %.3f',INdim+1,INdim,psi));
subplot(122);
    plot( x(:,1), x(:,2), '.' ); axis('image'); axis('off');
    title('First 2-dimensions only');
    disp('Points only reside on the unit circle if input data dimensionality is 1');

%% Surrogate testing

D = 1e3;        % surrogates
progress = 0;
psis = nan(D,1);
fprintf('\nSurrogate testing... %g surrogates of %g samples with dimensionality %g\n<',D,N,INdim+1);
for d = (1:D)
    currentprogress = round(100*d/D);
    if ( currentprogress > progress )
        progress = currentprogress;
        if (mod(progress,10)==0)
            if (progress~=100)
                fprintf('|');
            end;
        else
            fprintf('.');
        end;
    end;
    psis(d) = psin( 2*pi*rand(N,5) );
end;
fprintf('>\n');

pval = 1 - find( sort( psis ) >= psi, 1, 'first' )/length(psis);
if isempty(pval)
    pval = 0;
end;

figure;
    hist( psis ); hold('on');
    plot( psi*[1 1], ylim, 'r', 'linewidth', 2 );
    title(sprintf('PSI = %.3f, surrogate one-tailed p-value = %.3f',psi,pval));
