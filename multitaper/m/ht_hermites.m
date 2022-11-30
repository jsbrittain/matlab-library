function h = ht_hermites( N, K, rate, scale )

% Check input parameters
if (~exist('scale'))
    scale = [];
end;
% Defaults
if (isempty(scale))
    scale = 1;
end;

% Form hermite functions
t=(-(N/2-1):N/2)';              % Cover complete time axis
t=t/rate/scale;             % Time deformation to ellipical regions
H=zeros(length(t),K);
h=zeros(length(t),K);
for k = (0:K-1)
    % Form hermite polynomials (recursive structure)
    if (k==0),     H(:,1)=ones(length(t),1);
    elseif (k==1), H(:,2)=2*t;
    else           H(:,k+1)=2*t.*H(:,k)-2*(k-1)*H(:,k-1);
    end;
    % Generate hermite functions
    h(:,k+1)=H(:,k+1).*exp(-t.^2/2)/(pi^(1/4)*sqrt(2^k*gamma(k+1)));
end;
