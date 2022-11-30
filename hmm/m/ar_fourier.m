function sp = ar_fourier(a,freqs,rate,sigma2)

p=length(a);

sp = 1;
for k=(1:p)
    sp = sp + a(k)*exp(-1i*2*pi*freqs/(rate)*k);
end;
sp = sqrt(sigma2)./sp;
