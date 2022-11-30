function [spy11,spy22,dy,dyv]=jk_comp_sp(sp11,sp22,params);
%
% Function to compare spectra
%
% Applicable to the multi-sample case of unequal sample sizes.
% Use with jk_comp_psp.m for plotting.
%
% Method devised for ln but includes a switch to arbitrary log base,
% usually 10. (see note in Halliday 1995 for details).
%
% References
%   H.Bokil, D.J.Thomson & P.P.Mitra (2005), Comparing Spectra and Coherences
%       for Samples of Unequal Size.  Conference paper (689.22).
%

% Spectral comparison from jackknife estimates

% Determine corrective factor for alternative log bases
log_base=10;                % Usual for plotting spectra (use exp(1) for natural log)
const=1/log(log_base);      % equals 1 for natural log

% Determine DOF (2m,2n) for data (m,n = no. jackknife estimates)
m=size(params.jk11,2);
n=size(params.jk22,2);

% Calculate log spectra minus bias terms (all samples)
spy11=(log(sp11)-psi(m)+log(m))*const;                      % Digamma implemented by psi(X) equiv psi(0,X)
spy22=(log(sp22)-psi(n)+log(n))*const;

% Calculate log spectra minus bias terms (delete-one estimates)
y11=(log(params.jk11)-psi(m-1)+log(m-1))*const;
y22=(log(params.jk22)-psi(n-1)+log(n-1))*const;

% Determine dy (scaled diff of log spectra)
dy=(spy11-spy22)/sqrt(psi(1,m)*const^2+psi(1,n)*const^2);               % Trigamma implemented by psi(1,X)
dy1=(y11-spy22(:,ones(1,n)))/sqrt(psi(1,m-1)*const^2+psi(1,n)*const^2); % delete-one dy (ch.1)
dy1=m*dy(:,ones(1,m))-(m-1)*dy1;                                        %  => pseudovalues
dy2=(spy11(:,ones(1,m))-y22)/sqrt(psi(1,m)*const^2+psi(1,n-1)*const^2); % delete-one dy (ch.2)
dy2=n*dy(:,ones(1,n))-(n-1)*dy2;                                        %  => pseudovalues

% Calculate variance
meandy1=mean(dy1,2); meandy1=meandy1(:,ones(1,m));
meandy2=mean(dy2,2); meandy2=meandy2(:,ones(1,n));
dyv=mean((dy1-meandy1).^2,2)/(m-1)+mean((dy2-meandy2).^2,2)/(n-1);

%length(find((abs(dy))>1.96))/length(spy11)*100    % Measure of %age outside 95% (for rough calibration)
