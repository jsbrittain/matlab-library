function [coh,params]=jk_sp2_coh(sp11,sp22,sp12,params);
%
% Perform jack-knife error estimation for coherence
%
% Requires output from jack-knifing function mt_sp2.
%
% Output parameters
%   coh     Magnitude coherency (square to plot coherence)
%   params  Modified parameters structure (includes variance)
%

% atanh transform used for magnitude coherency, as (Bokil 2005),
%  but without the normalisation of (Thomson 1991).

% Jack-knife for coherence variance
M=size(params.jk12,2);      % Jackknife count
jkcoh=atanh(abs(params.jk12)./sqrt(params.jk11.*params.jk22));%-1/(2*params.jkcount-2);
meanjk=mean(jkcoh,2); meanjk=meanjk(:,ones(1,M));
params.jkcohv=(M-1)*mean((jkcoh-meanjk).^2,2);

% Estimate magnitude coherency
coh=atanh(abs(sp12)./sqrt(sp11.*sp22));%-1/(2*params.jkcount-2);
