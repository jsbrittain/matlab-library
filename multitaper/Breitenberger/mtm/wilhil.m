function [chi2]=wilhil(nu,p)
% Syntax: [chi2]=wilhil(nu,p);
% Wilson-Hilferty approximation to chi-squared.
% 
% Inputs: p - single-tailed probability
%         nu - degrees of freedom.
%
% The accuracy of the approximation is relatively
% poor for very small p and nu. The approximation
% is systematically low for small p and nu, and
% systematically large for large p and nu.
% 
% Accuracy will be within 1% for the following:
%     p = 0.005,  nu > 13. 
%     p = 0.01,   nu > 12. 
%     p = 0.025,  nu > 8. 
%     p = 0.05,   nu > 5.
%     p = 0.95,   nu > 1.
%
% More accurate results  (at the cost of speed)
% may be obtained using P. R. Shaw's 'chitable.m',
% available on the Mathworks ftp site.
%
% Either p or nu may be a vector, in which case 
% the output will be a matrix. 
% For example, wilhil(1:7, [.005 .05 .95 .995])
% will produce a table with 4 columns and 7 rows,
% containing the appropriate values of chi-squared.
%
% See e. g. Kendall and Stuart, Adv. Theory of 
% Statistics, Vol. 1, 1969.
%
% Written by Eric Breitenberger, version 10/8/95
% Please send comments and suggestions to eric@gi.alaska.edu

z=sqrt(2)*erfinv(2*p-1); % single-tailed
sig=2./(9*nu);
mu=1-sig;
sig=sqrt(sig);
if length(z)>1, 
  chi2=zeros(length(nu),length(z));
  for i=1:length(z)
    t1=nu.*(sig.*z(i)+mu).^3;
    chi2(:,i)=t1';
  end
else  
  chi2=nu.*(sig*z+mu).^3;
end
