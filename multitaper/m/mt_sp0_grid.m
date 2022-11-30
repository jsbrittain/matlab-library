function [sp11,sp22,sp12,params]=mt_sp0_grid(dat,expunge,duration,rate,opt_str)
%function [sp11,sp22,sp12,params]=mt_sp0_grid(dat,[expunge],duration,rate,opt_str)
%
% Generate spectra and coherence for grid plot
%
% Also accepts `expunge' vector (see mt_sp0_expunge.m)
%
%function [sp11,sp22,sp12,params]=mt_sp0_grid(dat,[expunge],duration,rate,opt_str)

% Determine input parameters
if (nargin<5)
    opt_str=rate;
    rate=duration;
    duration=expunge;
end;

% Traverse channels and generate spectra/bivariate parameters
M=size(dat,2);
for ch1=(1:M)
    for ch2=(ch1:M)
        disp([' Processing channel pair ' num2str(ch1) ' - ' num2str(ch2)]);
        if (nargin>4)
            [sp11{ch1,ch2},sp22{ch1,ch2},sp12{ch1,ch2},params{ch1,ch2}]=mt_sp0_expunge(dat(:,ch1),dat(:,ch2),expunge,duration,rate,opt_str);
        else
            [sp11{ch1,ch2},sp22{ch1,ch2},sp12{ch1,ch2},params{ch1,ch2}]=mt_sp0(dat(:,ch1),dat(:,ch2),duration,rate,opt_str);
        end;
    end;
end;
