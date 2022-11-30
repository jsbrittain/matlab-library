function [sp11,sp22,sp12,params]=mt_sp0_epochs_grid(epochs,rate,opt_str)
%function [sp11,sp22,sp12,params]=mt_sp0_epochs_grid(epochs,rate,opt_str)
%
% Generate spectra and coherence for grid plot
%
% Also accepts `expunge' vector (see mt_sp0_expunge.m)
%
%function [sp11,sp22,sp12,params]=mt_sp0_epochs_grid(epochs,rate,opt_str)

% Traverse channels and generate spectra/bivariate parameters
M=size(epochs,2);
for ch1=(1:M)
    for ch2=(ch1:M)
        disp([' Processing channel pair ' num2str(ch1) ' - ' num2str(ch2)]);
        [sp11{ch1,ch2},sp22{ch1,ch2},sp12{ch1,ch2},params{ch1,ch2}] = mt_sp0_epochs(epochs(:,ch1,:),epochs(:,ch2,:),rate,opt_str);
    end;
end;
