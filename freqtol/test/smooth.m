function out = smooth( in, width )

if ( exist('smooth','builtin') )
    out = builtin('smooth',in,width);
else
    out = runavg(in,width);
end
