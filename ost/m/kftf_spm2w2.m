function [sp,params]=kftf_spm2w2(dat1,dat2,trig,segoffset,segwidth,offset,duration,width,rate,seg_pwr,opt_str);
%function [sp,params]=kftf_spm2w2(dat1,dat2,trig,segoffset,segwidth,offset,duration,width,rate,seg_pwr,opt_str);
%
% Hybrid Kalman-Fourier analysis Time-Frequency (Type 2)
%
% Time-frequency analysis (perform analysis over offset sections of the
% trials).
%
% Used in pooled analysis routines - accepts vectors for offset/duration/width.
% Accepts a small segoffset vector for TFR style analyses and a segdur
% (0=start of offset location; given as portion of width).
%
%              |----------------------|
%            trig         ^________________ duration(ms)
%                  |---------------|
%               offset     ^_______________ width (both portions of duration)
%                  |-->|-->|-->|-->|
%              segoffset ^___^___^_________ segwidth (both portions of width; scalar or vector)
%
% Input parameters
%       dat1            Data vector channel 1
%       dat2            Data vector channel 2
%       trig            Trigger times
%       segoffset       Segment offset vector for TFR (portion of width)
%       segwidth        Segment duration vector for TFR (portion of width)
%       offset          Trial offset (portion of duration) (scalar or vector)
%       duration        Trial duration (msec) (scalar or vector)
%       width           Section width as a portion of duration (<=1)
%       rate            Sampling rate
%       seg_pwr         Segment power for Fourier transform
%       opt_str         Options string
%                           r<0|1|2>    Rectify channels
%                           W<W>        Multitaper analysis (using W; default 2.5)
%                           Q<q>        Process noise (default 0.1)
%                           f<fmax>     Maximum frequency (default Nyquist)
%                           s           Use Kalman smoother
%                           k<n>        KF method (0-null,1-kf,2-kflog(default),3-ff)
%                           j           Jackknife variance
%
% Output parameters
%       sp              Matrix of spectral components, cols: (3rd dim are trials)
%                           1           Auto-spectra ch.1
%                           2           Auto-spectra ch.2
%                           3           Cross-spectra ch.1,2
%                           4           Coherence
%       params          Parameters structure
%
% KF implementations
%       null            Return periodogram estimates (parameters structure will not fully correlate)
%       kf              Time-domain periodogram filtering (note distributional assumptions)
%       kflog           Log-domain (nonlinear) filtering
%       ff              Forgetting-factor (Q<%age> becomes delta(<1) value)
%
%function [sp,params]=kftf_spm2w2(dat1,dat2,trig,segoffset,segwidth,offset,duration,width,rate,seg_pwr,opt_str);

% Check input parameters
if (nargin~=11)
    error(' Incorrect number of input arguments');
end;
if (nargout~=2)
    error(' Incorrect number of output arguments');
end;
% if ((length(segwidth)>1) & (find(mean(segwidth)~=segwidth)) & (seg_pwr==0)) % Check widths not unequal for seg_pwr=0
%     error(' Segment widths cannot be of unequal portions with seg_pwr=0 (natural length).');
% end;
debug = true;

% Determine data parameters
sections=size(segoffset,2);

% Convert msecs to samples
duration=duration*rate/1000;
% Ensure possible scalar arguments converted to vectors
if (length(duration)~=length(trig))
    duration=duration*ones(length(trig),1);
end;
if (size(segoffset,1)~=length(trig))
    segoffset=segoffset(ones(length(trig),1),:);
end;
if (size(segwidth,1)~=length(trig))
    segwidth=segwidth(ones(length(trig),1),:);
end;
if (size(segwidth,2)~=sections)
    segwidth=segwidth(:,ones(1,sections));
end;

if (debug)
    disp(  'Variable sizes');
    disp([ '      offset (' num2str(size(offset)) ')']);
    disp([ '    duration (' num2str(size(duration)) ')']);
    disp([ '   segoffset (' num2str(size(segoffset)) ')']);
    disp([ '       width (' num2str(size(width)) ')']);
end;

% Perform KF-analysis within a time-slice
for ind=1:sections
    % Display progress
    disp(['Optimal spectral tracking - TFR (' int2str(ind) ' of ' int2str(sections) ')']);
    % Calculate offset and duration vectors (converted to msecs)
    offset2=round(offset.*duration + segoffset(:,ind).*width.*duration)*1000/rate;
    duration2=floor(segwidth(:,ind).*width.*duration)*1000/rate;
    % Perform analysis
    [sp{ind},params]=kf_spm2w(dat1,dat2,trig,offset2,duration2,rate,seg_pwr,[opt_str ' n']);
end;
if (length(sp)==1)
    sp=sp{:};
end;