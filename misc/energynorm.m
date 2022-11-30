function dat=energynorm(dat);

dat=dat/sqrt(sum(dat.^2));
