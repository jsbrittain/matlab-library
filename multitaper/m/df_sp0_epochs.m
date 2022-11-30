function [sp11,sp22,sp12,params,coh11,coh22,coh12]=df_sp0_epochs(dat1,dat2,rate,opt_str)
%[sp11,sp22,sp12,params,[coh11,coh22,coh12]]=df_sp0_epochs(dat1,dat2,rate,opt_str)
%
% Multitaper analysis
% Type 0
%
% Options as mt_sp2
%
% DO NOT USE PADDING OPTIONS
%
%[sp11,sp22,sp12,params,[coh11,coh22,coh12]]=df_sp0_epochs(dat1,dat2,rate,opt_str)

% Reconstruct segmented time series from epochs
dat1=squeeze(dat1);
dat2=squeeze(dat2);
duration=size(dat1,1)*1000/rate;
dat1=dat1(:);
dat2=dat2(:);

% Perform analysis over segments
[sp11,sp22,sp12,params,coh11,coh22,coh12]=df_sp0(dat1,dat2,duration,rate,opt_str);
