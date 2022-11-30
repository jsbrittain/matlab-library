% Test signal
rate = 100; K = 7.5*rate;
s = 0*sin(2*pi*5*(0:K-1)/rate + pi/2).*sin(2*pi*0.1*(0:K-1)/rate) + 0.1*randn(1,K);
omega0 = 2*pi*5;
M = 1;

mu0 = 0;%0.0001;      % Frequency adaptation (0=FLC)
mu1 = 0.01;         % FLC phase tracking parameter
muhat = 0.1;        % Amplitude tracking control

% Test FLC
dt = 1/rate;
[y,hy,omega] = wflc( s, omega0, M, mu0, mu1, muhat, 1/rate );

% Plot output
figure;
subplot(2,4,(1:3));
    plot(s); hold on; plot(y);
    plot( angle(hy)/10 );
    plot( abs(hy) );
    title(['M = ' num2str(M)]);
subplot(2,4,4);
    plot(omega/2/pi*rate);
    title('Frequency');
subplot(223);
    plot(s); hold on; plot(y,'linewidth', 2);
    plot( angle(hy)/10 );
    plot( abs(hy), 'linewidth', 2 );
    xlim([rate 3*rate]);
subplot(224);
    plot(s); hold on; plot(y,'linewidth', 2);
    plot( angle(hy)/10 );
    plot( abs(hy), 'linewidth', 2 );
    xlim([rate 1.5*rate]);

% Fwd-bwd
[yf,hyf,omegaf] = wflc( s, omega0, M, mu0, mu1, muhat, 1/rate );
yf = yf(:); hyf = hyf(:); omegaf = omegaf(:);
[yb,hyb,omegab] = wflc( flipud(s(:)), omega0, M, mu0, mu1, muhat, 1/rate );
yb = flipud(yb(:)); hyb = flipud(conj(hyb(:))); omegab = flipud(omegab(:));

figure;
subplot(3,1,(1:2));
    plot( s ); hold on;
    plot( yf ); plot( yb );
    plot( abs(hyf) ); plot( abs(hyb) );
    plot( (abs(hyf) + abs(hyb))/2 );
subplot(313);
    plot( angle(hyf) ); hold on; plot( angle(hyb) );
    plot( (angle(exp(1i*(angle(hyf)-angle(hyb))))), 'linewidth', 2 );
    plot( xlim, [0 0], 'k' )
