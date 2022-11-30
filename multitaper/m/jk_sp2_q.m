function [q,params]=jk_sp2_q(sp12,params);
%
% Perform jack-knife error estimation for cumulant density
%
% Requires output from jack-knifing function mt_sp2.
%
% Output parameters
%   q       Cumulant density
%   params  Modified parameters structure (includes variance)
%

% Determine if neurospec spectral norm used
q_norm=1;
dur=params.padcount*params.rate/1000;       % Duration samples
if (~isfield(params,'mtparams'))            % Required if called from mt_sp.m
    params.mtparams=params;
end;
% Neurospec normalisation
if (params.mtparams.spec_norm/(params.duration*params.rate/1000)==2*pi)
    q_norm=2*pi;
end;

% Jack-knife for cumulant variance (re-create full sp12 and ifft for cumulant)
M=size(params.jk12,2);                      % Jackknife count
jkq=q_norm*real(fftshift(ifft(params.jk12)));
meanjk=mean(jkq,2); meanjk=meanjk(:,ones(1,M));
params.jkqv=(M-1)*mean((jkq-meanjk).^2,2);

% Estimate cumulant
q=q_norm*real(fftshift(ifft(sp12)));
params.qlags=[-floor(dur/2):ceil(dur/2-1)]*1000/params.rate;
