function [permstats,channamesPerm,analysis] = analysisDriftPerm_SMRLoadTextMark( filenames, channames, chandesc, markerstrData, markerstrSham, frange, seglen, truncate, permutations )

overlapsecs = 0;
if (isempty(truncate))
    truncate = false;
end;

% Determine channel types
emgs = find( strcmp(chandesc,'EMG') );
not_emgs = setdiff( (1:length(chandesc)), emgs );

% Load data
preproc = [];
confrom_sampling_rates = true;
[hepochs,head] = SMRgetUsingMarkerData( filenames, markerstrData, channames, preproc, overlapsecs, confrom_sampling_rates );

% Pre-process settings -- EMGs
preprocEMG = [];
preprocEMG{1}.hpf.filtername = 'butter';
preprocEMG{1}.hpf.filterorder = 3;
preprocEMG{1}.hpf.cornerfreq = 50;
preprocEMG{1}.rectify = 1;

% Pre-process settings -- all channels including EMGs
preproc = [];
preproc{1}.lpf.filtername = 'butter';
preproc{1}.lpf.filterorder = 3;
preproc{1}.lpf.cornerfreq = frange(2);
preproc{1}.hpf.filtername = 'butter';
preproc{1}.hpf.filterorder = 3;
preproc{1}.hpf.cornerfreq = frange(1);
preproc{1}.hilbert = 1;

% Pre-process data
for k = (1:length(hepochs))
    % EMG preprocessing
    [hepochs{k}(:,emgs),head{k}(emgs)] = son_preproc( hepochs{k}(:,emgs), head{k}(emgs), preprocEMG );
    % All channel preprocessing
    [hepochs{k},head{k}] = son_preproc( hepochs{k}, head{k}, preproc );
end;

% Check that some data was loaded
if (isempty(hepochs))
    permstats = [];
    channamesPerm = [];
    analysis = [];
    return;
end;
% Check all records for dead channels
deadchans = false(1,length(head{1}));
for k = (1:length(head))
    deadchans = deadchans | [head{k}.dead];
end;
% Simplify header
head = head{1};

% Load sham data (after real data in-case there are no real markers in file)
[hepochsSham,headSham] = SMRgetUsingMarkerData( filenames, markerstrSham, channames, [], overlapsecs );
% Pre-process sham data
for k = (1:length(hepochsSham))
    [hepochsSham{k}(:,emgs),headSham{k}(emgs)] = son_preproc( hepochsSham{k}(:,emgs), headSham{k}(emgs), preprocEMG );
    [hepochsSham{k},headSham{k}] = son_preproc( hepochsSham{k}, headSham{k}, preproc );
end;
% Check all records for dead channels
for k = (1:length(headSham))
    deadchans = deadchans | [headSham{k}.dead];
end;
headSham = headSham{1};
deadchans( not_emgs ) = false;

% Cull dead channels from stim and sham datasets (do after load to simply channel orders)
for k = (1:length(hepochs)),     hepochs{k}     = hepochs{k}(:,~deadchans);     end;
for k = (1:length(hepochsSham)), hepochsSham{k} = hepochsSham{k}(:,~deadchans); end;
channames = channames(~deadchans);
chandesc = chandesc(~deadchans);
% Re-determine channel types
emgs = find( strcmp(chandesc,'EMG') );
not_emgs = setdiff( (1:length(chandesc)), emgs );

% Segment datasets
for k = (1:length(hepochs))
    hepochs{k} = segmentData( hepochs{k}, head, seglen );
end;
hepochs = cat(2,hepochs{:});
for k = (1:length(hepochsSham))
    hepochsSham{k} = segmentData( hepochsSham{k}, head, seglen );
end;
hepochsSham = cat(2,hepochsSham{:});

% Sort out channel order
channames = channames([ emgs not_emgs ]);
chandesc = chandesc([ emgs not_emgs ]);
ch = struct('accs',[],'emgs',[],'tacs',[]);
for k = (1:length(chandesc))
    ch.(lower(chandesc{k})) = find(strcmp(chandesc,chandesc{k}));
end;

% Replace sham tacs with frequency matched surrogate from stim data
freqtacs = median( angle(exp(1i*diff(angle( reshape(hepochs(:,:,ch.tacs),size(hepochs,1)*size(hepochs,2),1) )))) )*head(ch.tacs).rate/(2*pi);
if ( true )
    ph = 2*pi*freqtacs*(1:(size(hepochsSham,1)*size(hepochsSham,2)))/headSham(ch.tacs).rate;
    ph = reshape( ph, [ size(hepochsSham,1) size(hepochsSham,2) ] );
    hepochsSham(:,:,ch.tacs) = cos( ph ) + 1i*sin( ph );
end;

% Equalise lengths
if ( truncate )
    trials = min([ size(hepochs,2) size(hepochsSham,2) ]);
    hepochs = hepochs(:,1:trials,:);
    hepochsSham = hepochsSham(:,1:trials,:);
end;

% Call permutation analysis with segmented data
[permstats,channamesPerm,analysis] = analysisDriftPerm_segdata( hepochs, hepochsSham, ch, channames, permutations );
