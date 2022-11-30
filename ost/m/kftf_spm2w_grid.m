function [sp,params]=kftf_spm2w_grid(dat,trig,offset,duration,width,rate,seg_pwr,opt_str)
%function [sp,params]=kftf_spm2w_grid(dat,trig,offset,duration,width,rate,seg_pwr,opt_str)
%
% Generate spectra and coherence for grid plot
%
% Also accepts `expunge' vector (see mt_sp0_expunge.m)
%
%function [sp,params]=kftf_spm2w_grid(dat,trig,offset,duration,width,rate,seg_pwr,opt_str)

% Traverse channels and generate spectra/bivariate parameters
M=size(dat,2);
for ch1=(1:M)
    for ch2=(ch1:M)
        disp([' Processing channel pair ' num2str(ch1) ' - ' num2str(ch2)]);
        [sp{ch1,ch2},params{ch1,ch2}]=kftf_spm2w(dat(:,ch1),dat(:,ch2),trig,offset,duration,width,rate,seg_pwr,opt_str);
    end;
end;
