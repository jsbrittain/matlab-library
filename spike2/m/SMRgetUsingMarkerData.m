function [epoch,header] = SMRgetUsingMarkerData( filenames, markerstr, channames, preproc, overlapsecs, conform_sampling_rates )

if (~exist('conform_sampling_rates','var'))
    conform_sampling_rates = [];
end;
if (isempty(overlapsecs))
    overlapsecs = 0;
end;
if (~iscell(preproc))
    preproc = { preproc };
end;
if (isempty(conform_sampling_rates))
    conform_sampling_rates = false;
end;

% Parse files and report markers
epoch = cell(0); header = cell(0);
for fileno = (1:length(filenames))
    fid = fopen( filenames{fileno}, 'r' );
    chanlist = SONChanList(fid);
    chanix = find([chanlist.kind]==8);
    if (isempty(chanix))
        fclose(fid);
        continue;
    end;
    chan = chanlist(chanix).number;
    
    % Isolate markers
    [data,head] = SONGetTextMarkerChannel(fid, chan);
    fclose(fid);
    data.markers = int8( data.markers - double('0') );
    data.text = mat2cell( data.text, ones(1,size(data.text,1)), size(data.text,2) );
    data.text = cellfun( @(x) deblank(char(x)), data.text, 'UniformOutput', false );
    if (~all(strcmp(data.text(1:2:end),data.text(2:2:end))))
        error(['Unmatched Text Markers: ' filenames{fileno}]);
    end;
    
    % Check if required data is in file
    kk = false(1,length(data.text));
    for k = (1:length(kk))
        kk(k) = any( strcmpi(data.text{k},markerstr) );
    end;
    if ~any( kk )
        continue;
    else
        % Remove unwanted markers
        data.timings = data.timings(kk,:);
        data.markers = data.markers(kk,:);
        data.text    = data.text(kk,:);
    end;
    
    % Read data
    [datk,headk] = son_getdata( filenames{fileno}, channames, preproc );
    if (isempty(datk))
        continue;
    end;
    time = (1:size(datk,1))/headk(1).rate;
    
    % Recurse markers
    for segk = (1:length(data.text)/2)
        kk = round((data.timings(2*segk-1)-overlapsecs)*headk(1).rate):round((data.timings(2*segk)+overlapsecs)*headk(1).rate);
        kk(kk<1) = []; kk(kk>size(datk,1)) = [];        %%% Add NaNs instead %%%
        epoch{end+1} = datk( kk, : );
        header{end+1} = headk;
    end;
end;

% Conform sampling rates between files
if ( conform_sampling_rates && (length(header)>1))
    if any(diff(cellfun(@(x) x(1).rate,header)))
        % Conform sampling rate --- need to read data first to know resampling frequency
        warning('Inconsistent sampling rates detected... resampling...');
        preproc{end+1}.resample.rate = header{1}(1).rate;
        [epoch,header] = SMRgetUsingMarkerData( filenames, markerstr, channames, preproc, overlapsecs );
    end;
end;
