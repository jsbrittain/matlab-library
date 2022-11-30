function trig = trigThreshold( x, level )

trig = find( diff( x > level )>0 );
