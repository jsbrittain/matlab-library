function dat=powernorm(dat);

dat=dat/sqrt(sum(dat.^2)/length(dat));
