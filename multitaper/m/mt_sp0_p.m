function [sp11,sp22,sp12,params]=mt_sp0_p(dat1,dat2,duration,rate,opt_str)
%function [sp11,sp22,sp12,params]=mt_sp0_p(dat1,dat2,duration,rate,opt_str);
%
% Multi-taper analysis type 0
%
% Implementation:
%   One time-series and one point process (spike train)
%
% Input parameters
%   Data parameters
%       dat1        Spike times ch.1
%       dat2        Time-series ch.2
%       trig        Trigger
%       offset      Offset (ms)
%       duration    Duration (ms)
%       rate        Sampling rate
%   Analysis parameters
%       frange      Frequency range (may be [] for defaults)
%       opt_str     Options string (as mt_sp2.m)
%
%function [sp11,sp22,sp12,params]=mt_sp0_p(dat1,dat2,duration,rate,opt_str);

% Convert times of occurance to a 0-1 process
N=size(dat2,1);
sp_dat1=zeros(N,1);                 % Assign variable space
sp_dat1(dat1)=1;                    % Populate with spikes

if (max(dat1)>N)
    error(' Spike times exceed length of time-series.');
end;

% Analyse using time series routines (Mean subtracted trials in mtm routine)
[sp11,sp22,sp12,params]=mt_sp0(sp_dat1,dat2,duration,rate,opt_str);
