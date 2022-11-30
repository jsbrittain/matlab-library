function [f,t,cl,M]=kf_fourierw(dat1,dat2,rate,seg_pwr,opt_str,S);
%function [f,t,cl,M]=kf_fourierw(dat1,dat2,rate,seg_pwr,opt_str,S);
%
% Segmented Fourier analysis type 0
% Moving window segmentation
%
% For use with KF analysis routines
%
% First five parameters as sp2a2_m.
% Additional parameters
%   S   Number of segments per group
%
%function [f,t,cl,M]=kf_fourierw(dat1,dat2,rate,seg_pwr,opt_str,S);

% Check input parameters
if (isempty(opt_str))
    opt_str=' ';    % Compatibility check for sp2a2_ma
end;

% Determine data parameters
N=length(dat1);

% Determine group and segment parameters
T=2^seg_pwr;                % Segment length (T)
L=fix(N/T);                 % Number of complete segments (L)
R=L*T;                      % Total number of samples (R=LT)
G=S*T;                      % Group length (G)
M=L-S+1;                    % Number of (complete) groups

% Perform segmented analysis
for ind=1:M
    % Segment data into (possibly incomplete) groups
    tt=(T*(ind-1)+1):(T*(S+ind-1));
    segdat1=dat1(tt);
    segdat2=dat2(tt);
    
    % Perform spectral analysis
    disp(['Group ' int2str(ind) ' of ' int2str(M)]);
    [f(:,:,ind),t(:,:,ind),cl(ind)]=sp2a2_ma(segdat1,segdat2,rate,seg_pwr,opt_str);
end;
