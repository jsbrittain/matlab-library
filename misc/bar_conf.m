function bar_conf( x, e )

bar( x - e );
hold on;
for n=1:length(x)
    plot( n*[1 1], x(n) + e(n)*[-1 +1], 'k' );
    plot( n+0.05*[-1 1], x(n) + e(n)*[1 1], 'k' );
end;
