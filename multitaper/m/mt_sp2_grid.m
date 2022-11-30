function [sp11,sp22,sp12,params]=mt_sp2_grid(dat,trig,offset,duration,rate,opt_str)
%function [sp11,sp22,sp12,params]=mt_sp2_grid(dat,trig,offset,duration,rate,opt_str)
%
% Generate spectra and coherence for grid plot
%
%function [sp11,sp22,sp12,params]=mt_sp2_grid(dat,trig,offset,duration,rate,opt_str)

% Traverse channels and generate spectra/bivariate parameters
M=size(dat,2);
for ch1=(1:M)
    for ch2=(ch1:M)
        disp([' Processing channel pair ' num2str(ch1) ' - ' num2str(ch2)]);
        [sp11{ch1,ch2},sp22{ch1,ch2},sp12{ch1,ch2},params{ch1,ch2}]=mt_sp2(dat(:,ch1),dat(:,ch2),trig,offset,duration,rate,opt_str);
    end
end
