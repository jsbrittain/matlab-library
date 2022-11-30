function [dat,time,rate]=vitaportascii_load(filename)

chcount=6;

% Load channels
for ind=(1:chcount)
    fid=fopen([filename '.a0' int2str(ind)],'r');
    
    % header lines
    for ind2=(1:10), fgetl(fid); end;
    
    % Load channel data
    dat{ind}=fscanf(fid,'%g',[Inf]);
    
    fclose(fid);
end;

time=[];
rate=[];
