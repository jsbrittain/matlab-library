function erpplot_image(erptime,epochs,wind,ix)
%function erpplot_image(erptime,epochs,wind)
%
%
%
%function erpplot_image(erptime,epochs,wind)

% ERP Plot

ep=squeeze(epochs);
out=zeros(size(ep,2),size(ep,1));

% sort according to `ix'
if (exist('ix'))
    [sortedix,ix]=sort(ix);
    ep = ep(:,ix);
end;

% Smooth
for n=(1:size(ep,1))
    out(:,n)=smooth(ep(n,:),wind);
end;

figure;
subplot(6,1,[1:5]);
imagesc(erptime,(wind:size(out,1)-wind+1),out(wind:end-wind+1,:));
subplot(6,1,6);
erpplot(erptime,ep);
