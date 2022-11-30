function [f,t,cl,out1,out2]=kf_fourierw(dat1,dat2,rate,seg_pwr,opt_str,S);
%function [f,t,cl,M]=kf_fourierw(dat1,dat2,rate,seg_pwr,opt_str,S);
%function [f,t,cl,sp,M]=kf_fourierw(dat1,dat2,rate,seg_pwr,opt_str,S);
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
%function [f,t,cl,[sp,]M]=kf_fourierw(dat1,dat2,rate,seg_pwr,opt_str,S);

% Check output requirements
sp_out=0;
if (nargout>4)
    sp_out=1;
end;

% Perform options in opt_str once on entire data segment
[dat1,dat2,opt_str]=kf_dataopt(dat1,dat2,opt_str);

% Determine data parameters
N=length(dat1);

% Determine group and segment parameters
T=2^seg_pwr;                % Segment length (T)
L=fix(N/T);                 % Number of complete segments (L)
R=L*T;                      % Total number of samples (R=LT)
G=S*T;                      % Group length (G)
M=L-S+1;                    % Number of (complete) groups

% Calculate maximum matrix requirements
seg_size=2^seg_pwr;
f_max=seg_size/2+1;
t_max=seg_size;
sp_max=seg_size/2;
f=zeros(f_max,5,M);
t=zeros(t_max,2,M);
sp=zeros(sp_max,4,M);

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

% Remove duplicate frequency component at end of f matrix
% (only occurs with sp2a2_ma)
f=f(1:(size(f,1)-1),:,:);

% Output spectral coefficients if required
if (sp_out)
    out1=sp;
    out2=M;
else
    out1=M;
end;
