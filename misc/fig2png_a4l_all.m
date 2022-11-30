% Convert all fig files in directory to equivalently named jpg files

contents = dir('*.fig');
imageformat = 'png';

for ind=1:size(contents,1)
    fh=open(contents(ind).name);
    namelen = size(contents(ind).name,2);
    jpgname=contents(ind).name;
    jpgname(namelen-2:namelen)='png';
    fig2image_a4l(fh, jpgname, imageformat);
    close(fh);
end;
