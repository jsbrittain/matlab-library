function out = low_pass(data, threshold);

length = size(data,1);

for n=1:length
    if (data(n)>threshold)
        out(n)=0;
    else
        out(n)=data(n);
    end
end

out=out';