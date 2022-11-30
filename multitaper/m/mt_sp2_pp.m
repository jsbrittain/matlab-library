function [sp11,sp22,sp12,params]=mt_sp2_pp(dat1,dat2,trig,offset,duration,rate,opt_str);
%function [sp11,sp22,sp12,params]=mt_sp2_pp(dat1,dat2,trig,offset,duration,rate,opt_str);
%
% Multi-taper analysis type 2
%
% Implementation:
%   Point process (spike train) data
%
% Input parameters
%   Data parameters
%       dat1        Spike times ch.1
%       dat2        Spike times ch.2
%       trig        Trigger
%       offset      Offset (ms)
%       duration    Duration (ms)
%       rate        Sampling rate
%   Analysis parameters
%       frange      Frequency range (may be [] for defaults)
%       opt_str     Options string
%                       r<0,1,2>    rectify (none, ch1, ch1&2)
%                       W           multitaper bandwidth (default: 5)
%                       f           maximum analysis frequency (default: Nyquist)
%                       t           trial-averaging only (no multitaper)
%                       j<m>        jackknife estimates
%                                     (m=method,0:trials,1:tapers,default:0)
%
%function [sp11,sp22,sp12,params]=mt_sp2_pp(dat1,dat2,trig,offset,duration,rate,opt_str);

% Convert times of occurance to a 0-1 process
N=max([max(dat1) max(dat2)]);
sp_dat1=zeros(N,1);                 % Assign variable space
sp_dat2=zeros(N,1);                 %  |
sp_dat1(dat1)=1;                    % Populate with spikes
sp_dat2(dat2)=1;                    %  |

% Analyse using time series routines (Mean subtracted trials in mtm routine)
[sp11,sp22,sp12,params]=mt_sp2(sp_dat1,sp_dat2,trig,offset,duration,rate,opt_str);
