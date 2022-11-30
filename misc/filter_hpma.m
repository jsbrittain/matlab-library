function out = filter_hpma(dataset, points)
% function out = filter_hpma(dataset, points)
%
% High Pass Moving Average Filter
%
% Inputs
%    dataset    Source process
%    points     Number of points to average over
%

if mod(points, 2) ~= 1
    error(' Input dataset must be of the form 2n+1');
end

length = size(dataset,1);
n = (points-1)/2;

for i=(n+1):(length - n)
    diff(i-n) = mean(dataset((i-n):(i+n)));
end

out = dataset((n+1):(length-n)) - diff';
