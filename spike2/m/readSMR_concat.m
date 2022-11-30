function [hdat,head,ch,pcacoeff,latent] = readSMR_concat( filename, channames, chandesc, frange, varargin )
%
% Input parameters
%   filename    Cell array of filenames to concatenate
%   channames   Cell array of channel names to load
%   frange      Band-pass range
%
% Parameter value pairs
%   pcacoeff    Use specified PCA rotation matrix
%

% Input parser
p = inputParser;
addRequired(p,'filename',@(x) iscell(x) || ischar(x));
addRequired(p,'channames',@iscell);
addRequired(p,'chandesc',@iscell);
addRequired(p,'frange',@(x) (isnumeric(x)) && ((length(x)==2) || (isempty(x))));
addParameter(p,'pcacoeff',[],@isnumeric);
addParameter(p,'rate',[],@(x) ((isnumeric(x) && isscalar(x)) || (isempty(x))));
addParameter(p,'isAmpMod',false);
parse( p, filename, channames, chandesc, frange, varargin{:} );

% Recover parameter-value inputs from parser
pcacoeff = p.Results.pcacoeff;
newrate = p.Results.rate;
isAmpMod = p.Results.isAmpMod;
latent = [];

% Form cell array if filename is a string
if (~iscell(filename))
    filename = {filename};
end;

% Concatenate records
for fileno = (1:length(filename))
    % Get data
    preproc = [];
    [datk,headk] = son_getdata( filename{fileno}, channames, preproc );
    if (size(datk,2)~=length(channames))
        error('Data length does not match provided channel names!');
    end;
    ch = [];
    ch.accs = find(strcmpi(chandesc,'Acc'));
    ch.emgs = find(strcmpi(chandesc,'EMG'));
    ch.tacs = find(strcmpi(chandesc,'TACS') | strcmpi(chandesc,'TACSOUT'));
    % Add PCA'd accelerometers
    ch.accs_pca = length(headk)+(1:length(ch.accs));
    if (~isempty(ch.accs_pca))
        % High-pass accelerometer at lower tremor frequency before PCA'ing
        [bh,ah] = butter(3,frange(1)/(headk(ch.accs(1)).rate/2),'high');
        accf = nan(size(datk,1),length(ch.accs_pca));
        for k = (1:length(ch.accs_pca))
            accf(:,k) = filtfilt(bh,ah,datk(:,ch.accs(k)));
            headk(ch.accs_pca(k)) = headk(ch.accs(k));
            headk(ch.accs(k)).title = sprintf('PCA-%g',k);
        end;
        if (isempty(pcacoeff))
            % Use data-centric rotation
            [pcacoeff,datk(:,ch.accs_pca),latent] = pca( accf );
        else
            % Use specified rotation (already high-pass filtered so mean subtracted)
            datk(:,ch.accs_pca) = (pcacoeff*accf')';
        end;
        clear('accf');
    end;
    % Preprocess (rectify) EMG
    preproc = [];
    preproc.hpf.filterorder = 3;
    preproc.hpf.filtername = 'butter';
    preproc.hpf.cornerfreq = 55;
    preproc.abs = [];
    [datk(:,ch.emgs),headk(ch.emgs)] = son_preproc( datk(:,ch.emgs), headk(ch.emgs), preproc );
    % High-pass and rectify Amplitude-Modulated stimulation waveform
    if ( isAmpMod )
        disp('Treating stim as amTACS...');
        [bh1,ah1] = butter( 3, 1.0/(headk(ch.tacs).rate/2), 'high' );
        abs( filtfilt(bh,ah,datk(:,ch.tacs)) );
    end;
    % Bandpass filter all data
    preproc = [];
    preproc.lpf.filterorder = 3;
    preproc.lpf.filtername = 'butter';
    preproc.lpf.cornerfreq = frange(2);
    preproc.hpf.filterorder = 3;
    preproc.hpf.filtername = 'butter';
    preproc.hpf.cornerfreq = frange(1);
    preproc.dcremove = [];
    % Hilbert transform all data
    preproc.hilbert = [];
    [hdatk,headk] = son_preproc( datk, headk, preproc );
    hdatk = hdatk(headk(1).rate:(end-headk(1).rate),:);
    clear('datk');
    % Downsample (already bandpass filtered)
    if (~isempty(newrate))
        hdatk0 = hdatk; hdatk = []; dt = 1/newrate;
        for k = (1:size(hdatk0,2))
            hdatk(:,k) = interp1( (1:size(hdatk0,1))/headk(k).rate, hdatk0(:,k), (dt:dt:(size(hdatk0,1)/headk(k).rate)) );
            headk(k).rate = newrate;
        end;
    end;
    % Concatenate with previous datasets
    if (fileno==1)
        hdat = hdatk;
        head = headk;
    else
        hdat = [ hdat; hdatk ];
        head = [ head headk ];
    end;
    clear('datk','headk');
end;
