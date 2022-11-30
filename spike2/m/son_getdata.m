function [dat,head] = son_getdata( filename, channames, preproc )
% Assumes channels are uniform length, if not, load using separate calls

if (~exist('preproc','var'))
    preproc = struct([]);
end;
if (~iscell(preproc))
    preproc = {preproc};
end;

allow_empty_channels = false;
if (isfield(preproc{1},'allow_empty_channels'))
    allow_empty_channels = preproc{1}.allow_empty_channels;
end;

% Form cell array from string (scalar)
if (~iscell(channames))
    channames = {channames};
end;

% Deconstruct multiple entries
if any(cellfun(@iscell,channames,'UniformOutput',true))
    preproc1 = [];
    preproc1.allow_empty_channels = true;
    for k = (1:length(channames))
        [datk,headk] = son_getdata( filename, channames{k}, preproc1 );
        if (~isempty(datk))
            dat(:,k) = datk;
            head(k) = headk;
        end;
    end;
    [dat,head] = son_preproc( dat, head, preproc );
    return;
end;

% Open file
fid = fopen(filename,'r');
if (fid==-1)
    error('Invalid filenames!');
end;

% Find channel
casesensitive = false;
try
    chans = cellfun( @(x) SONChanFindTitle( fid, x, false ), channames, 'UniformOutput', true );
catch
    if ( (~allow_empty_channels) && any(cellfun(@isempty,cellfun( @(x) SONChanFindTitle( fid, x, false ), channames, 'UniformOutput', false ))))
        warning('Unusual situation - Cannot find requested channels --- skipping file!');
        dat = [];
        head = [];
        return;
    else
        chankk = find( cellfun( @(x) ~isempty( SONChanFindTitle( fid, x, false ) ), channames, 'UniformOutput', true ) );
        chans = cellfun( @(x) SONChanFindTitle( fid, x, false ), channames(chankk), 'UniformOutput', true );
    end;
end;

% Load channel data
chanlist = SONChanList( fid );
switch (chanlist(find([chanlist.number]==chans(1),1,'first')).kind)
    case 1, [dat,head] = SONGetADCChannel( fid, chans(1), 'scale' );                % Normal ADC
    case 9, [dat,head] = SONGetChannel( fid, chans(1), 'scale' );                   % Alek spliced
    otherwise, error(' Unexpected channel type!');
end;
for v = (2:length(chans))
    % Read data
    switch (chanlist(find([chanlist.number]==chans(v),1,'first')).kind)
        case 1, [dat0,head(v)] = SONGetADCChannel( fid, chans(v), 'scale' );        % Normal ADC
        case 9, [dat0,head(v)] = SONGetChannel( fid, chans(v), 'scale' );           % Alek spliced
        otherwise, error(' Unexpected channel type!');
    end;
    dat(:,v) = [ dat0(1:min(size(dat,1),length(dat0))); dat0(end)*ones(max(0,length(dat)-size(dat0,1)),1) ];
end;
dat = double(dat);

% Convert sampling interval to sampling rate
for v = (1:length(chans))
    switch (chanlist(find([chanlist.number]==chans(v),1,'first')).kind)
        case 1, head(v).rate = round(1e6/head(v).sampleinterval);
        case 9, head(v).rate = round(1e6/head(v).sampleinterval);
        otherwise, error(' Unexpected channel type!');
    end;
end;

% Close file
fclose(fid);

% Mark dead channels
for v = (1:length(chans))
    head(v).dead = false;
    % Select dead channel criterion
    switch ( 2 )
        case 1,     % Look for data at (scaled) channel maximum --- may need to check validity of this!
            criterion = ( ( quantile(dat(:,v)/head(v).scale,0.025) <= -2.5 ) && ( quantile(dat(:,v)/head(v).scale,0.975) >= +2.5 ) );
        case 2,     % Look for a percentage of data at its maximum (or minimum) value
            clip_limit = 0.10;      % Percentage of clipped data to permit
            criterion = ( mean( (dat(:,v)==max(dat(:,v))) | (dat(:,v)==min(dat(:,v))) ) >= clip_limit );
        otherwise
            error('Unknown dead channel selection criterion selected!');
    end;   
    if ( criterion )
        head(v).dead = true;
    end;
end;

% Preprocessing
[dat,head] = son_preproc( dat, head, preproc );
