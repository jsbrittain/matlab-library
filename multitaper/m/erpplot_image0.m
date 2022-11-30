function erpplot_image0(erptime,epochs,wind)
%function erpplot_image0(erptime,epochs,wind)
%
%
%
%function erpplot_image0(erptime,epochs,wind)

% ERP Plot

ep=squeeze(epochs);
out=zeros(size(ep,2),size(ep,1));
for n=(1:size(ep,1))
    out(:,n)=smooth(ep(n,:),wind);
end;

imagesc(erptime,(wind:size(out,1)-wind+1),out(wind:end-wind+1,:));
