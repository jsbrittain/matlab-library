% Test signal
rate = 100; K = 10*rate; clear s;
%s = sin(2*pi*5*(0:K-1)/rate + pi/4).*sin(2*pi*0.1*(0:K-1)/rate) + 0.1*randn(1,K);
omega0 = 2*pi*5/rate;
omega = omega0 + rand(1)*2*pi;
for k = 1:K
    omega0 = omega0 + 0.000*randn(1);
    freq(k) = omega0;
    omega = omega + omega0;
    s(k) = sin(omega) + 0.01*randn(1);
end;
s = s.*sin(2*pi*0.1*(0:K-1)/rate);
omega0 = 2*pi*5;
M = 1;

mu = 0.01;

% Test FLC
frange = [4 6]; G = 10;
[y,hy] = blflc( s, frange, G, mu, 1/rate );

% Plot output
subplot(2,4,(1:3));
    plot(s); hold on; plot(y);
    title(['M = ' num2str(M)]);
    
    plot( abs( hy ) );
    plot( angle( hy )/10 );

yf = blflc( s(:), frange, G, mu, 1/rate );
yf = yf(:);
yb = blflc( flipud(s(:)), frange, G, mu, 1/rate );
yb = flipud(yb(:));

subplot(4,1,3);
    plot( s ); hold on;
    plot( yf );
    plot( yb );
%subplot(4,1,4);
%    plot( angle() ); hold on;
%    plot( angle() );
