function out = filter_rec(dataset)
% function out = filter_rec(dataset)
%
% Simple Recursive filter
%

len = length(dataset);
out = dataset;  % Create structure before computation
                % - speeds computation exponentially
for i=2:len
    out(i,1) = dataset(i) - dataset(i-1);
end
