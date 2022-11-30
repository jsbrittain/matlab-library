function [sp11,sp22,sp12,params]=mt_sp0_expunge(dat1,dat2,expunge,duration,rate,opt_str)
%[sp11,sp22,sp12,params]=mt_sp0_expunge(dat1,dat2,expunge,duration,rate,opt_str)
%
% Multitaper analysis
%   Type 0
%
% Expunge routine permtits the additional vector `expunge' which rejects
% all samples marked `0' (false) and prevents boundary effects by selecting
% epochs within the accepted range only.
%
% Options as mt_sp2.m
%
%[sp11,sp22,sp12,params]=mt_sp0_expunge(dat1,dat2,expunge,duration,rate,opt_str)

% Determine segments
offset=0;                       % Default offset
dt=duration*rate/1000;          % msecs -> samples

% Determine triggers                                        (CRUDE METHOD FOR QUICK IMPLEMENTATION)
trig=(1:dt:length(dat1)-dt);    % Disjoint sections
reject=find(~expunge);
ind=1;
while (ind<=length(trig))       % Reject epochs in expunged zone
    if (~isempty(find(reject>trig(ind) & reject<(trig(ind)+dt),1)))
        trig=[trig(1:(ind-1)) trig((ind+1):end)];
    else
        ind=ind+1;
    end;
end;

% Perform analysis over segments
[sp11,sp22,sp12,params]=mt_sp2(dat1,dat2,trig,offset,duration,rate,opt_str);
