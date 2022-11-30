function [sp11,sp22,sp12,params]=mt_sp0_pp(dat1,dat2,duration,rate,opt_str)
%[sp11,sp22,sp12,params]=mt_sp0_pp(dat1,dat2,duration,rate,opt_str);
%
% Multitaper analysis
% Type 0
%
% Options as mt_sp2
%
%[sp11,sp22,sp12,params]=mt_sp0_pp(dat1,dat2,duration,rate,opt_str);

% Determine segments
N=max([max(dat1) max(dat2)]);
offset=0;                       % Default offset
dt=duration*rate/1000;          % msecs -> samples
trig=[1:dt:N-dt+1];             % Disjoint sections

% Perform analysis over segments
[sp11,sp22,sp12,params]=mt_sp2_pp(dat1,dat2,trig,offset,duration,rate,opt_str);
