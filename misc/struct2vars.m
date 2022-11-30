function struct2vars(s)

fnames = fieldnames(s);
fvals = struct2cell(s);

for i=1:length(fnames)
    assignin('caller',fnames{i},fvals{i});
end
