function [sp11,sp22,sp12,params]=mt_sp2_welch(dat1,dat2,trig,offset,duration,rate,seglen,overlap,opt_str)
%function [sp11,sp22,sp12,params]=mt_sp2_welch(dat1,dat2,trig,offset,duration,rate,seglen,overlap,opt_str)
%
% Wrapper function for mt_sp2 implementing basic windowing within the
% region of interest.  Will eventually perform Welch's overlapping
% segments, but for now just segments into disjoint sections and performs
% multitaper analysis on each.
%
% Additional input parameters
%       seglen          Segment length (msecs)
%       overlap         Percentage overlap [0 100)
%
%function [sp11,sp22,sp12,params]=mt_sp2_welch(dat1,dat2,trig,offset,duration,rate,seglen,overlap,opt_str)

if (overlap~=0)
    error(' Function not implemented for overlapping segments yet.');
end;

% Reduce passed duration argument and segment region by adding stepped
% trigger points (DISJOINT ONLY AT PRESENT)
trig2=trig; trig=[];
for ind=(1:length(trig2))
    trig=[trig; trig2(ind)+(1:seglen*rate/1000:duration*rate/1000)];       % DOES NOT ACCOUNT FOR INCOMPLETE SEGMENTS
end;
duration=seglen;

% Perform type 2 multitaper analysis
[sp11,sp22,sp12,params]=mt_sp2(dat1,dat2,trig,offset,duration,rate,opt_str);
