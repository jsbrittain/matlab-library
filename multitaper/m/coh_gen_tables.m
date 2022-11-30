%% Generate coherence tables with specified precision

reps = 100;
N = 1e5;                    % Data length (double this for FFT then take +ve freqs only)
Mlist = (2:2:20);             % List of trials to generate coherence for

alphalist = (0:0.01:1);     % Coherence
cohvals = (0:0.01:1);       % Evaluate at these coherence values

% Convert desired coherence (alpha) to alpha-parameter for mixing process
for n = (1:length(alphalist))
    c = alphalist(n);
    r = roots([ c 0 2*c-4 0 c ]);
    alphalist(n) = real( r(end) );
end;

% Reserve memory space
coh = zeros(length(cohvals),length(alphalist),length(Mlist));

for rep = (1:reps)
    for k = (1:length(alphalist))

        % Update display
        if (k>1), toc, end;
        disp([ 'Rep ' num2str(rep) ' / ' num2str(reps) ', Alpha ' num2str(k) ' / ' num2str(length(alphalist)) ]);
        tic

        % For each coherence value
        for n = (1:length(Mlist))
            M = Mlist(n);

            % Bivariate time-series
            x0 = randn(2*N,M);
            y0 = randn(2*N,M);

            % Mixing process (see Lacheaux)
            x = x0 + alphalist(k)*y0;
            y = alphalist(k)*x0 + y0;

            % Take FFT
            fx = fft(x,[],1); fx = fx(1:end/2,:);
            fy = fft(y,[],1); fy = fy(1:end/2,:);

            % Construct spectra
            Sxy = fx.*conj(fy);
            Sxx = abs(fx).^2;
            Syy = abs(fy).^2;

            % Histogram of coherence over frequencies (which are independent)
            coh(:,k,n) = coh(:,k,n) + hist( ( abs(mean(Sxy,2)).^2./(mean(Sxx,2).*mean(Syy,2)) ), cohvals ).'/size(Sxx,1)/reps;
            %coh(:,k,n) = coh(:,k,n)/sum(coh(:,k,n));

        end;
    end;
    coh0 = coh;
end;

clear('Sxy','Sxx','Syy','fx','fy','x','y','x0','y0','M','n','k','c','r');
save('cohlookup');

%% Non-parametric response - Test

if ( false )
    estimated = 0.2;

    % For a given K, determine most likely `estimated' coherence value
    m = 1;
    plot( alphalist, coh( dsearchn( cohvals', estimated ), :, m ) );

    [cohmag, cohpos] = max(coh( dsearchn( cohvals', estimated ), :, m ));
    title( [ 'Given estimated coherence = ' num2str(estimated) ' for K = ' num2str(Mlist(m)) ', likely real coherence = ' num2str(cohvals(cohpos)) ' (p = ' num2str(cohmag) ')' ] );
end;
