function cl=cp_conflim(X,p)

B=10000;
Y=zeros(B,1);
for ind=(1:B)
    Y(ind)=mean(X(round(rand(length(X)-1,1)+1)));
end;

[pdf,x]=hist(Y,length(X));
cdf=cumsum(pdf);
cdf=cdf/cdf(end);

cl=x(dsearchn(cdf',[p (1-p)]'));
