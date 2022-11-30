function fx=mmae_dens(z, mu, C);
%function fx=normal_d(z, mu, C);
%
% Adapted Normal (Gaussian) distribution
%
% z=observation value
% mu='a priori' state prediction
% C=estimate of error residuals
%

fx=inv(sqrt(2*pi*abs(C))) * exp( -0.5 * (z-mu)'*inv(C)*(z-mu) );
