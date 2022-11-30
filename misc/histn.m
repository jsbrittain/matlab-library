function [count,bins] = histn( x, bins )

if (~exist('bins'))
    bins = 20;
end;

[count,bins] = hist( x, bins );
count = count/sum(count(:));

if (nargout==0)
    try
        bar( bins, count );
    catch
        [count,bins] = hist( x, bins );
    end;
end;
