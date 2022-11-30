function ch = SONChanFindUnits( fid, str, casesensitive )

if (~exist('casesensitive'))
    casesensitive = [];
end;
if (isempty(casesensitive))
    casesensitive = true;
end;

chanlist = SONChanList( fid );
for k = (1:length(chanlist))
    chanlist(k).units = getfield(SONChannelInfo(fid,chanlist(k).number),'units');
end;
if (casesensitive)
    ch = find( strcmp({chanlist.units},str), 1, 'first' );
else
    ch = find( strcmpi({chanlist.units},str), 1, 'first' );
end;

if ~isempty( ch )
    ch = chanlist(ch).number;
end;
