function out = mexican_hat(t)

out = (1 - t.^2) .* exp(-1 * t.^2 / 2);
