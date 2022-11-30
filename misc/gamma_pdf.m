function pdf=gamma_dist(x,k,theta);

pdf=x.^(k-1).*exp(-x/theta)/theta^k/gamma(k);
