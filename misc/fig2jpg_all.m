% Convert all fig files in directory to equivalently named jpg files

contents = dir('*.fig');

for ind=1:size(contents,1)
    fh=open(contents(ind).name);
    namelen = size(contents(ind).name,2);
    jpgname=contents(ind).name;
    jpgname(namelen-2:namelen)='jpg';
    fig2jpeg_a4l(fh, jpgname);
    close(fh);
end;
