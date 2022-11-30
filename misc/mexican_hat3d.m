function out = mexican_hat3d(x, y)

out = (2-(x.^2+y.^2)).*exp(-(x.^2+y.^2)/2);
