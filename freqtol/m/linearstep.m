function y = linearstep( coeffs, x )
%
% coeffs = [ m c ]
%

m = -coeffs(1);
c =  coeffs(2);

%  y  =  m.(x - c)
y = m.*(x - c);
