function pdf=lognormal(x,mu,sigma);

pdf=exp(-(log(x-mu)).^2/(2*sigma^2))./(x*sigma*sqrt(2*pi));
