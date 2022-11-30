function [sp,params] = kftf_spm2w2_epochs(dat1,dat2,segoffset,segwidth,width,rate,seg_pwr,opt_str)

% Remove singleton dimensions
dat1=squeeze(dat1);
dat2=squeeze(dat2);

% Extract epoch parameters
offset = 0;                         % Offset (null)
dur=size(dat1,1);                   % Duration of epoch (samples)
duration = dur/rate*1000;           %  | (msecs)
M = size(dat1,2);                   % Trial count

% Flatten time-series
dat1=dat1(:);
dat2=dat2(:);

% Determine trigger locations
trig = (1:dur:M*dur)';

% Perform kftf_spm2w analysis
[sp,params] = kftf_spm2w2(dat1,dat2,trig,segoffset,segwidth,offset,duration,width,rate,seg_pwr,opt_str);
