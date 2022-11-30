function [f,t,cl,out1,out2]=kf_fourier2(dat1,dat2,trig,offset,duration,rate,seg_pwr,opt_str,S);
%function [f,t,cl,M]=kf_fourier2(dat1,dat2,trig,offset,duration,rate,seg_pwr,opt_str,S);
%function [f,t,cl,sp,M]=kf_fourier2(dat1,dat2,trig,offset,duration,rate,seg_pwr,opt_str,S);
%
% Segmented Fourier analysis type 2
%
% For use with KF analysis routines
%
% First five parameters as sp2a2_m.
% Additional parameters
%   S   Number of segments per group
%
%function [f,t,cl,[sp,]M]=kf_fourier2(dat1,dat2,trig,offset,duration,rate,seg_pwr,opt_str,S);

% Perform options in opt_str once on entire data segment
%[dat1,dat2,opt_str]=kf_dataopt(dat1,dat2,opt_str);

% Check output requirements
sp_out=0;
if (nargout>3)
    sp_out=1;
end;

% Determine data parameters
N=length(dat1);
offset=offset*rate/1000;
duration=duration*rate/1000;

% Determine group and segment parameters
L=length(trig);
M=ceil(L/S);

% Calculate maximum matrix requirements
seg_size=2^seg_pwr;
f_max=seg_size/2+1;
t_max=seg_size;
sp_max=seg_size/2;
f=zeros(f_max,5,M);
t=zeros(t_max,2,M);
sp=zeros(sp_max,4,M);

% Perform spectral analysis on the M groups of S trials
for ind=1:M
    % Extract triggered data to pass to analysis routines
    % (more efficient than passing it all)
    newtrig=trig((S*(ind-1)+1):min(S*ind,L));
    segdat1=dat1((newtrig(1)+offset):(newtrig(length(newtrig))+offset+duration));
    segdat2=dat2((newtrig(1)+offset):(newtrig(length(newtrig))+offset+duration));
    newtrig=newtrig-newtrig(1)-offset;
    
    % Perform analysis over group
    disp(['Segment group ' int2str(ind) ' of ' int2str(M)]);
    [f(:,:,ind),t(:,:,ind),cl(ind),sp(:,:,ind)]=sp2a2_ma(segdat1,segdat2,newtrig,offset,duration,rate,seg_pwr,opt_str);
end;

% Remove duplicate frequency component at end of f matrix
% (only occurs with sp2a2_ma)
f=f(1:(size(f,1)-1),:,:);

% Combine third and fourth spectral columns into a complex cross-spectra
sp(:,3,:)=sp(:,3,:)+i*sp(:,4,:);
sp=sp(:,1:3,:);

% Output spectral coefficients if required
if (sp_out)
    out1=sp;
    out2=M;
else
    out1=M;
end;
