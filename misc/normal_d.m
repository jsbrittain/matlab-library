function fx=normal_d(x, mu, sigma);
%function fx=normal_d(x, mu, sigma);
%
% Normal (Gaussian) distribution
%
% x=distribution point
% mu=mean
% sigma=variance ^(1/2)
%

fx=(1/(sqrt(2*pi)*sigma)) * exp( -((x-mu).^2)/(2*sigma^2) );
