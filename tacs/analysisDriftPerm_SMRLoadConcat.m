function [permstats,channamesPerm,analysis,pcacoeff,latent] = analysisDriftPerm_SMRLoadConcat( filename, filenameSham, channames, chandesc, frange, seglen, dsrate, truncate, permutations, use_same_rotation, bincount, isAmpMod, sham_mode, splitmode )

if (~exist('bincount','var'))
    bincount = [];
end;
if (~exist('isAmpMod','var'))
    isAmpMod = [];
end;
if (~exist('sham_mode','var'))
    sham_mode = [];
end;
if (~exist('splitmode','var'))
    splitmode = [];
end;

if (isempty(isAmpMod))
    isAmpMod = false;
end;
if (isempty(truncate))
    truncate = false;
end;
if (isempty(sham_mode))
    sham_mode = 1;          % 1 = Use TACS frequency from stim in sham
end;
if (isempty(splitmode))
    splitmode = 0;
end;

% Load data
[hdatSham,headSham,ch,pcacoeff{1},latent{1}] = readSMR_concat( filenameSham, channames, chandesc, frange, 'rate', dsrate, 'isAmpMod', isAmpMod );
if (use_same_rotation), pcacoeff{2} = pcacoeff{1}; else pcacoeff{2} = []; end;
[hdat,head,~,pcacoeff{2},latent{2}] = readSMR_concat( filename, channames, chandesc, frange, 'pcacoeff', pcacoeff{2}, 'rate', dsrate, 'isAmpMod', isAmpMod );
channames = {head.title};

% Implement legacy pre-processing
legacy  = struct([]);
fnames = fieldnames(legacy);
if (~isempty(legacy))
    warning(' !!! LEGACY options are ACTIVE !!!');
    disp(fnames);
    for k = (1:length(fnames))
        switch (fnames{k})
            case {'flip'},     % Flip accelerometer
                hdat = abs(hdat).*exp(1i*( pi+angle(hdat) ));      % Flip acc
                warning(' --- Flipping accelerometer channels...');
            case {'runnorm'},     % 30 secs running amplitude normalisation
                hdat(:,3) = (abs(hdat(:,3))./runavg(abs(hdat(:,3)),30*head(1).rate)).*exp(1i*angle(hdat(:,3)));
                warning(' --- Running average power normalisation [30 secs]...');
        end;
    end;
end;

% Replace sham tacs with frequency matched surrogate from stim data
if ( true )
    switch ( sham_mode )
        case 0,     % Leave alone (good for comparing different stim sessions)
        case 1,     % Use stim frequency for sham
            freqtacs = median( angle(exp(1i*diff(angle(hdat(:,ch.tacs))))) )*head(ch.tacs).rate/(2*pi);
            ph = 2*pi*freqtacs*(1:size(hdatSham,1))/headSham(ch.tacs).rate;
            hdatSham(:,ch.tacs) = cos( ph ) + 1i*sin( ph );
        case 2,     % Use *adjusted* (stim-sham) relative frequency for sham
            freqtacs = median( angle(exp(1i*diff(angle(hdat(:,1))))) )*head(1).rate/(2*pi);
            warning('using *adjusted* tacs frequency for sham!');
            ph = 2*pi*freqtacs*(1:size(hdatSham,1))/headSham(ch.tacs).rate;
            hdatSham(:,ch.tacs) = cos( ph ) + 1i*sin( ph );
        otherwise
            error('Unknown sham mode.');
    end;
end;

% Segment datasets
hdat = segmentData( hdat, head, seglen );
hdatSham = segmentData( hdatSham, headSham, seglen );

% Subset selection of epochs
switch ( splitmode )
    case 0, % Do nothing
    case 1, % Low median split
        ampSeg = nanmean( nanmean( abs( hdat ), 1), 3 );
        hdat = hdat(:,ampSeg<median(ampSeg),:);
        ampSegSham = nanmean( nanmean( abs( hdatSham ), 1), 3 );
        hdatSham = hdatSham(:,ampSegSham<median(ampSegSham),:);
    case 2, % High median split
        ampSeg = nanmean( nanmean( abs( hdat ), 1), 3 );
        hdat = hdat(:,ampSeg>=median(ampSeg),:);
        ampSegSham = nanmean( nanmean( abs( hdatSham ), 1), 3 );
        hdatSham = hdatSham(:,ampSegSham>=median(ampSegSham),:);
end;

% Equalise lengths
if ( truncate )
    trials = min([ size(hdat,2) size(hdatSham,2) ]);
    hdat = hdat(:,1:trials,:);
    hdatSham = hdatSham(:,1:trials,:);
end;

% Call permutation analysis with segmented data
[permstats,channamesPerm,analysis] = analysisDriftPerm_segdata( hdat, hdatSham, ch, channames, permutations, bincount );
