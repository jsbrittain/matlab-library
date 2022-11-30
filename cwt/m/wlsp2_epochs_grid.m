function [Wsp,wlparam]=wlsp2_epochs_grid(epochs,rate,df,frange,opt_str,mother,wlopt)
%function [Wsp,wlparam]=wlsp2_epochs_grid(epochs,rate,df,frange,opt_str,mother,wlopt)


for ch1 = (1:size(epochs,2))
    for ch2 = (ch1:size(epochs,2))
        disp([ ' Pair ' num2str(ch1) ' - ' num2str(ch2) ', (of ' num2str(size(epochs,2)) ')' ]);
        [Wsp{ch1,ch2},wlparam{ch1,ch2}] = wlsp2_epochs(epochs(:,ch1,:),epochs(:,ch2,:),rate,df,frange,opt_str,mother,wlopt);
    end;
end;
