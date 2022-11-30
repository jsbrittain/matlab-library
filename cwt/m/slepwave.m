function [psi,lambda]=slepwave(M,nwin,p,pc);
%function [psi,lambda]=slepwave(M,nwin,p,pc);
%
% Slepian wavelets
%
% Input parameters
%       M           Data points (default: 100)
%       nwin        Mumber of real-valued wavelets to generate (default: 6)
%       p           Time-bandwidth product (default: 2.5)
%       pc          Time-bancentre product (default: 5.0)
%
% Output parameters
%       psi         Slepian wavelets
%       lambda      Energy concentrations within the time-freq area
%
% Generation code lifted from appendix of Lily & Park 1995.
%
% Ref: Lilly & Park (1995) "Multiwavelet spectral and polarization analyses
%       of seismic records". Geophys.J.Int. 122, pp.1001-1021.
%
%function [psi,lambda]=slepwave(M,nwin,p,pc);

% Assign default parameters if unspecified
if (nargin<4)
    pc=5.0;
end;
if (nargin<3)
    p=2.5;
end;
if (nargin<2)
    nwin=6;
end;
if (nargin<1)
    M=100;
end;

% Interpolate if M large (>Minterp)
Minterp=200;
if (M>Minterp)
    % Interpolate Minterp point vectors to M points
    [psi,lambda]=slepwave(Minterp,nwin,p,pc);
    dM=Minterp/M;
    psi=interp1(0:Minterp-1,psi,0:dM:Minterp-dM,'spline');
else
	% Calculate Toeplitz matrix row
	fw=p/M;
	fo=pc/M;
	x=linspace(2*pi,(M-1)*2*pi,M-1);
	sink1=2.0*sin((fo+fw)*x)./x;
	sink1=[2.0*(fo+fw) sink1];
	sink2=2.0*sin((fo-fw)*x)./x;
	sink2=[2.0*(fo-fw) sink2];
	sink=sink1-sink2;
	slep=toeplitz(sink);
	[V,D]=eig(slep);
	
	% Order eigenvalues
	[eigv,k]=sort(diag(D));
	V=V(:,k);
	
	% Reorder last nwin eigenvectors and eigenvalues to first lambda=eigv(n:-1:n-nwin+1)
	lambda=eigv(M:-1:M+1-nwin);
	psi=V(:,M:-1:M+1-nwin);
	
	% Set sign conventions for wavelets
	for k=1:nwin
        if (psi(2,k)<0)
            psi(:,k)=-psi(:,k);
        end;
	end;
end;

% Clear lambda if not output
if (nargout<2)
    clear('lambda');
end;
