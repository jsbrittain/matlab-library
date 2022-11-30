function [dat,time,head]=smr_load(filename,trange,chanlist)
%function [dat,time,head]=smr_load(filename)
%
% Load SMR channels
%
%function [dat,time,head]=smr_load(filename)

warning off
addpath ~/matlab/library/thirdparty/sigtool/sigTOOL/s'igTOOL Neuroscience Toolkit'/File/menu_Import/'group_NeuroScience File Formats'/son
addpath ~/matlab/library/thirdparty/sigtool/sigTOOL/s'igTOOL Neuroscience Toolkit'/File/menu_Import/'group_NeuroScience File Formats'/son/SON32
warning on

disp(' Loading Spike2 file ...');

if (~exist('trange'))
    trange=[];
end;

% Open file
fid=fopen(filename);
header=SONFileHeader(fid);

% Load channels
chans=SONChanList(fid);
if (exist('chanlist'))
    chans=chans(chanlist);
end;
disp(['  ' num2str(length(chans)) ' channels']);
for ind=(1:length(chans))
    [dat{ind},head{ind}]=SONGetChannel(fid,chans(ind).number);
    time{ind}=(head{ind}.start:head{ind}.sampleinterval*1e-6:head{ind}.stop);
    
    if (~isempty(trange))
        pts=dsearchn(time{ind}',trange');
        dat{ind}=dat{ind}(pts(1):pts(2));
        time{ind}=time{ind}(pts(1):pts(2));
    end;
end;

% Close file
fclose(fid);
disp(' done.');
