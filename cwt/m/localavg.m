function y=localavg(x,m);

N=length(x);
r=(m-1)/2;
y=zeros(size(x));
y(r+1)=sum(x(1:m));
for t=r+1:N-r-1
    y(t+1)=y(t)-x(t-r)+x(t+r+1);
end;
