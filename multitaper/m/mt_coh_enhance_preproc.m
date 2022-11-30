function unew = mt_coh_enhance_preproc( u, v, a, b, rate )

if (~exist('a','var'))
    a = [];
end;
if (~exist('b','var'))
    b = [];
end;

if (isempty(a))
    a = hamming(rate+1);
end;
if (isempty(b))
    b = hamming(8*rate+1);
end;

u = u(:);
v = v(:);
a = a(:);
b = b(:);

K = (length(a)-1)/2;
M = (length(b)-1)/2;

unew = zeros(size(u));
for k = (-K:K)
    [ k 2*K+1 ]
    
    uv = u.*circshift(v,k);
    cuvb = conv( uv, b, 'same' );
    for n = ((K+1):(length(u)-K))
        unew(n) = unew(n) + a(K+k+1)* cuvb(n) *v(n-k);
    end;
end;
