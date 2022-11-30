function y = trilinearstep( coeffs, x )
%
% coeffs = [ m1 m2 c1 c2 ]
%

m1 = -coeffs(1);
m2 = -coeffs(2);
c1 =  coeffs(3);
c2 =  coeffs(4);

%  y  =  /   m1(x-c1)        f <= c1
%        |   m2(x-c2)        f >  c2
%        \   0               otherwise
y = (x<=c1).*(m1*(x-c1)) + (x>c2).*(m2*(x-c2));
