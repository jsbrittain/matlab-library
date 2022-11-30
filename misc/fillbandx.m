function fillbandx(x1,x2,h,colstr,alpha)

if (~exist('h'))
    h=gca;
end;
if (~exist('colstr'))
    colstr='g';
end;
if (~exist('alpha','var'))
    alpha=0.25;
end;

hold(h,'on');
y1=min(ylim(h)); y2=max(ylim(h));

for n=(1:length(x1))
    fill([x1(n) x2(n) x2(n) x1(n) x1(n)],[y2 y2 y1 y1 y2],colstr,'edgecolor','none','facealpha',alpha);
end;
