function [ph,params]=jk_sp2_ph(sp12,params);
%
% Perform jack-knife error estimation for phase
%
% Requires output from jack-knifing function mt_sp2.
%
% Output parameters
%   ph      Phase
%   params  Modified parameters structure (includes variance)
%

% Jack-knife for cumulant variance (re-create full sp12 and ifft for cumulant)
M=size(params.jk12,2);                      % Jackknife count
jkph=angle(params.jk12);
meanjk=mean(jkph,2); meanjk=meanjk(:,ones(1,M));
params.jkphv=(M-1)*mean((jkph-meanjk).^2,2);

% Estimate cumulant
ph=angle(sp12);
