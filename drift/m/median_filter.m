function y = median_filter( x, n, display_progress, skip )

if (~exist('display_progress'))
    display_progress = true;
end;
if (~exist('skip'))
    skip = 1;
end;
progress = 0;

fprintf('\n Median filter: ');
y = median(x)*ones(size(x,1)/skip,1);
for k=(n:skip:(length(x)-n))
    % Display progress
    if ( display_progress )
        if (round(100*k/length(x))>progress)
            progress = round(100*k/length(x));
            if (mod(progress,10)==0)
                fprintf('%g',progress);
            else
                fprintf('.');
            end;
        end;
    end;
    % Median filter
    tt=k+(-n/2:n/2);
    y(k+(0:(skip-1)))=median( x(tt) );
end;
fprintf('100%%\n');

if ( skip > 0 )
    switch ( 2 )
        case 1,
            y = smooth( y, 2*skip );
        case 2,
            t = n:skip:(length(x)-n);
            y = interp1( t, y(t), 1:length(x), 'spline' );
    end;
end;
