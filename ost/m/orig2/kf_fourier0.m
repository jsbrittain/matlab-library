function [f,t,cl,sp,M]=kf_fourier0(dat1,dat2,rate,seg_pwr,opt_str,S);
%function [f,t,cl,sp,M]=kf_fourier0(dat1,dat2,rate,seg_pwr,opt_str,S);
%
% Segmented Fourier analysis type 0
% Disjoint segregation
%
% For use with KF analysis routines
%
% First five parameters as sp2a2_m.
% Additional parameters
%   S   Number of segments per group
%
%function [f,t,cl,sp,M]=kf_fourier0(dat1,dat2,rate,seg_pwr,opt_str,S);

% Determine data parameters
N=length(dat1);

% Determine group and segment parameters
T=2^seg_pwr;                % Segment length (T)
L=fix(N/T);                 % Number of complete segments (L)
R=L*T;                      % Total number of samples (R=LT)
group_length=S*T;           % Group length
M=ceil(L/S);                % Number of (possibly incomplete) groups

% Perform segmented analysis
for ind=1:M
    % Segment data into (possibly incomplete) groups
    tt=(group_length*(ind-1)+1):min((group_length*ind),R);
    segdat1=dat1(tt);
    segdat2=dat2(tt);
    
    % Perform spectral analysis
    disp(['Segment group ' int2str(ind) ' of ' int2str(M)]);
    [f(:,:,ind),t(:,:,ind),cl(ind),sp(:,:,ind)]=sp2a2_ma(segdat1,segdat2,rate,seg_pwr,opt_str);
end;

% Remove duplicate frequency component at end of f matrix
% (only occurs with sp2a2_ma)
f=f(1:(size(f,1)-1),:,:);

% Combine third and fourth spectral columns into a complex cross-spectra
sp(:,3,:)=sp(:,3,:)+i*sp(:,4,:);
sp=sp(:,1:3,:);
