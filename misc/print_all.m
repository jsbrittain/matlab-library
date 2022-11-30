% Print all fig files in directory

contents = dir('*.fig');

for ind=1:size(contents,1)
    fh=open(contents(ind).name);
    print(fh);
    close(fh);
end;
