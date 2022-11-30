function [sp,M]=kf_fouriersp2w(dat1,dat2,trig,offset,duration,rate,seg_pwr,opt_str,S);
%function [sp,M]=kf_fouriersp2w(dat1,dat2,trig,offset,duration,rate,seg_pwr,opt_str,S);
%
% Segmented Fourier analysis type 2
% Moving window segmentation
%
% For use with KF analysis routines
%
% First five parameters as sp2a2_m.
% Additional parameters
%   S   Number of segments per group
%
%function [sp,M]=kf_fouriersp2w(dat1,dat2,trig,offset,duration,rate,seg_pwr,opt_str,S);

% Perform options in opt_str once on entire data segment
[dat1,dat2,opt_str]=kf_dataopt(dat1,dat2,opt_str);

% Determine data parameters
N=length(dat1);
offset=offset*rate/1000;
duration=duration*rate/1000;

% Determine group and segment parameters
L=length(trig);
M=L-2;

% Perform spectral analysis on the M groups of S trials
for ind=1:M
    % Extract triggered data to pass to analysis routines
    % (more efficient than passing it all)
    newtrig=trig(ind:(ind+S-1));
    segdat1=dat1((newtrig(1)+offset):(newtrig(length(newtrig))+offset+duration));
    segdat2=dat2((newtrig(1)+offset):(newtrig(length(newtrig))+offset+duration));
    newtrig=newtrig-newtrig(1)-offset;
    
    % Perform analysis over group
    disp(['Segment group ' int2str(ind) ' of ' int2str(M)]);
    [f_dummy,t_dummy,cl_dummy,sp(:,:,ind)]=sp2a2_ma(segdat1,segdat2,newtrig,offset,duration,rate,seg_pwr,opt_str);
end;
