% Test signal
rate = 100; K = 10*rate; clear s;
%s = sin(2*pi*5*(0:K-1)/rate + pi/4).*sin(2*pi*0.1*(0:K-1)/rate) + 0.1*randn(1,K);
omega0 = 2*pi*5/rate;
omega = omega0 + 10*rand(1)*2*pi;
for k = 1:K
    omega0 = omega0 + 0.000*randn(1);
    freq(k) = omega0;
    omega = omega + omega0;
    s(k) = sin(omega) + 0.01*randn(1);
end;
s = s.*sin(2*pi*0.1*(0:K-1)/rate)
s = s + 0.1*randn(size(s));
omega0 = 2*pi*5;
M = 5;

mu = 0.1;

% Test FLC
[y,hy] = flc( s, omega0, M, mu, 1/rate );

% Plot output
subplot(2,4,(1:3));
    plot(s); hold on; plot(y);
    plot( angle(hy)/10 );
    plot( abs(hy) );
    title(['M = ' num2str(M)]);
    %plot(abs(hilbert(y)));
subplot(2,4,4);
    plot( freq/2/pi*rate );
subplot(223);
    plot(s); hold on; plot(y,'linewidth', 2);
    plot( angle(hy)/10 );
    plot( abs(hy), 'linewidth', 2 );
    %xlim([1000 3000]);
subplot(224);
    plot(s); hold on; plot(y,'linewidth', 2);
    plot( angle(hy)/10 );
    plot( abs(hy), 'linewidth', 2 );
    %xlim([1000 1500]);
