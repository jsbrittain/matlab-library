function sp = armaspectrum(ar,ma,freqs,rate,sigma2)

p=length(ar);
q=length(ma);

spar = 1;
for k=(1:p)
    spar = spar + ar(k)*exp(-1i*2*pi*freqs/rate*k);
end;
spma = 1;
for k=(1:q)
    spma = spma + ma(k)*exp(-1i*2*pi*freqs/rate*k);
end;

sp = sigma2.*(abs(spma).^2)./(abs(spar).^2);
