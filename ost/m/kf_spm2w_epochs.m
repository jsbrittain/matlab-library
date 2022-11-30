function [sp,params] = kf_spm2w_epochs(dat1,dat2,rate,seg_pwr,opt_str)
%function [sp,params] = kf_spm2w_epochs(dat1,dat2,rate,seg_pwr,opt_str)
%
%function [sp,params] = kf_spm2w_epochs(dat1,dat2,rate,seg_pwr,opt_str)

% Remove singleton dimensions
dat1=squeeze(dat1);
dat2=squeeze(dat2);

% Extract epoch parameters
offset = 0;                         % Offset (null)
dur = size(dat1,1);                 % Duration of epoch (samples)
duration = dur/rate*1000;           %  | (msecs)
M = size(dat1,2);                   % Trial count

% Flatten time-series
dat1=dat1(:);
dat2=dat2(:);

% Determine trigger locations
trig = (1:dur:M*dur);

% Perform kf_spm2w analysis
[sp,params] = kf_spm2w(dat1,dat2,trig,offset,duration,rate,seg_pwr,opt_str);
