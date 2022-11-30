function out = filter_exp_ma(dataset, alpha, recursive)
% function out = filter_exp_ma(dataset, alpha)
%
% Exponentially Weighted Moving Average Filter
% (First Order Low Pass Filter)
%
% Input:
%   dataset     Time series
%   alpha       Filter constant alpha=n/(n+1) implying (1-alpha) = 1/(n+1)
%   recursive   Number of back steps to feed over
%

% fx_k = alpha * fx_k-1 + (1-alpha) * x_k

% To calculate value out(k):
out(k) = dataset(k)*(1-alpha) + dataset(k-1)*alpha*(1-alpha) + dataset(k-2)*alpha^2*(1-alpha) + dataset(k-3)*alpha^3;
