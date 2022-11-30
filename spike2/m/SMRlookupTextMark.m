function [data,head] = SMRlookupTextMark( filenames )

% Parse files and report markers
for fileno = (1:length(filenames))
    fid = fopen( filenames{fileno}, 'r' );
    chanlist = SONChanList(fid);
    chanix = find([chanlist.kind]==8);
    if (isempty(chanix))
        continue;
    end;
    chan = chanlist(chanix).number;
    
    [data(fileno),head(fileno)] = SONGetTextMarkerChannel(fid, chan);
    data(fileno).markers = int8( data(fileno).markers - double('0') );
    data(fileno).text = mat2cell( data(fileno).text, ones(1,size(data(fileno).text,1)), size(data(fileno).text,2) );
    data(fileno).text = cellfun( @(x) deblank(char(x)), data(fileno).text, 'UniformOutput', false );
    fclose(fid);
end;

if (~exist('data','var'))
    data = [];
    head = [];
end;
