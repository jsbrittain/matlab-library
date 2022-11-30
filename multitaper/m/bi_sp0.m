function [p11,bs11,p22,bs22,p12,bs12,bs21,params,bicoh11,bicoh22,bicoh12,bicoh21]=bi_sp0(dat1,dat2,duration,rate,opt_str)
%[p11,bs11,p22,bs22,p12,bs12,bs21,params,[bicoh11,bicoh22,bicoh12,bicoh21]]=bi_sp0(dat1,dat2,duration,rate,opt_str)
%
% Multitaper analysis
% Type 0
%
% Options as mt_sp2
%
%[p11,bs11,p22,bs22,p12,bs12,bs21,params,[bicoh11,bicoh22,bicoh12,bicoh21]]=bi_sp0(dat1,dat2,duration,rate,opt_str)

% Determine segments
offset=0;                       % Default offset
dt=duration*rate/1000;          % msecs -> samples
trig=[1:dt:length(dat1)-dt+1];  % Disjoint sections

% Perform analysis over segments
[p11,bs11,p22,bs22,p12,bs12,bs21,params,bicoh11,bicoh22,bicoh12,bicoh21]=bi_sp2(dat1,dat2,trig,offset,duration,rate,opt_str);
