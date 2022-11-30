function pdf=normpdf(x,mu,Sigma);
%function pdf=normpdf(x,mu,Sigma);
%
% Multivariate Normal Probability Distribution Function
%
% Input parameters
%       x       Matrix of multivariate random vectors
%                   (N x T) where N is multivariate dimensionality
%                             and T is the number of multivariate vectors
%       mu      Column vector of multivariate mean
%       Sigma   Covariance matrix (N x N)
%
%function pdf=normpdf(x,mu,Sigma);

N=length(mu);
pdf=zeros(N,size(x,2));
for ind=1:length(x)
    pdf(:,ind)=(1/((2*pi)^(N/2))*sqrt(det(Sigma)))*exp(-.5*transpose(x(:,ind)-mu)*inv(Sigma)*(x(:,ind)-mu));
end;

% MULTIVARIATE ABOVE DOES NOT SEEM TO WORK PROPERLY
pdf=(1/sqrt(2*pi*Sigma))*exp(-(x-mu).^2./(2*Sigma));
