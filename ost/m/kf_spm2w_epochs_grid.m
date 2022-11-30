function [sp,params] = kf_spm2w_epochs_grid(dat,rate,seg_pwr,opt_str)
%function [sp,params] = kf_spm2w_grid(dat,rate,seg_pwr,opt_str)
%
% Generate spectra and coherence for OST grid plot
%
%function [sp,params] = kftf_spm2w_grid(dat,rate,seg_pwr,opt_str)

% Extract epoch parameters
offset = 0;                         % Offset (null)
dur = size(dat,1);                  % Duration of epoch (samples)
duration = dur/rate*1000;           %  | (msecs)

% Determine trigger locations
Ntrials = size(dat,3);
trig = (1:dur:Ntrials*dur);

% Traverse channels and generate spectra/bivariate parameters
M=size(dat,2);
for ch1=(1:M)
    for ch2=(ch1:M)
        disp([' Processing channel pair ' num2str(ch1) ' - ' num2str(ch2)]);
        [sp{ch1,ch2},params{ch1,ch2}] = kf_spm2w( ...
            reshape(dat(:,ch1,:),size(dat,1)*size(dat,3),1),...
            reshape(dat(:,ch2,:),size(dat,1)*size(dat,3),1),...
            trig,offset,duration,rate,seg_pwr,opt_str);
    end;
end;
