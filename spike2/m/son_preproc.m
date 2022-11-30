function [dat,head] = son_preproc( dat, head, preproc )
% Assumes channels are uniform length, if not, load using separate calls

if (iscell(dat))
    for ind = (1:length(dat))
        [dat{ind},head{ind}] = son_preproc( dat{ind}, head{ind}, preproc );
    end;
end;
if (length(size(dat))>2)
    error('Data dimensionality too high - expecting vectors of channels');
end;

if (isempty(preproc))
    return;
end;
if (~iscell(preproc))
    preproc = {preproc};
end;
if all(cellfun(@isempty,preproc))
    return;
end;
chans = (1:size(dat,2));

% Preprocessing - Acc
for ind = (1:length(preproc))
    fnames = fieldnames(preproc{ind});
    for k = (1:length(fnames))
        switch ( fnames{k} )
            
            case 'lpf',         % Low-pass filter
                if any(diff([head.rate])), perchannelfilters = true;
                else                       perchannelfilters = false; end;
                for v = (1:length(chans))
                    if (perchannelfilters || (v==1))
                        [b,a] = feval(preproc{ind}.(fnames{k}).filtername,preproc{ind}.(fnames{k}).filterorder,preproc{ind}.(fnames{k}).cornerfreq/(head(v).rate/2),'low');
                    end;
                    dat(:,v) = filtfilt( b, a, filtfilt( b, a, dat(:,v) ) );
                end;
                
            case 'hpf',         % High-pass filter
                if any(diff([head.rate])), perchannelfilters = true;
                else                       perchannelfilters = false; end;
                for v = (1:length(chans))
                    if (perchannelfilters || (v==1))
                        [b,a] = feval(preproc{ind}.(fnames{k}).filtername,preproc{ind}.(fnames{k}).filterorder,preproc{ind}.(fnames{k}).cornerfreq/(head(v).rate/2),'high');
                    end;
                    dat(:,v) = filtfilt( b, a, filtfilt( b, a, dat(:,v) ) );
                end;
                
            case {'abs','rectify'},         % Absolute (full-wave rectification)
                dat = abs( dat );
                
            case 'hilbert',     % Hilbert transform
                for v = (1:length(chans))
                    dat(:,v) = hilbert( dat(:,v) );
                end;
                
            case {'dcremove','demean'},     % DC remove
                for v = (1:length(chans))
                    dat(:,v) = dat(:,v) - nanmean(dat(:,v));
                end;
                
            case {'detrend'},               % Detrend
                for v = (1:length(chans))
                    dat(:,v) = dat(:,v) - polyval( polyfit(1:size(dat,1),dat(:,v),1), 1:size(dat,1) );
                end;
                
            case {'resample'},              % Resample
                newrate = preproc{ind}.(fnames{k}).rate;
                dt = 1/newrate;
                if any([head.rate]~=newrate)
                    dat0 = dat; dat = [];
                    for v = (1:length(chans))
                        if ( newrate == head(v).rate )
                            dat(:,v) = dat0(:,v);
                            continue;
                        end;
                        if ( newrate < head(v).rate )
                            disp('   Applying anti-aliasing filter during resampling...');
                            [bl,al] = butter(3,(newrate/2.5)/(head(v).rate/2),'low');
                            dat0(:,v) = filtfilt( bl, al, dat0(:,v) );
                        end;
                        dat(:,v) = interp1( (1:size(dat0,1))/head(v).rate, dat0(:,v), (dt:dt:(size(dat0,1)/head(v).rate)), 'linear' );
                        head(v).rate = newrate;
                        head(v).interval = 1/newrate;
                    end;
                end;
                
            case {'fftnotch'},              % FFT notch
                for v = (1:length(chans))
                    dat(:,v) = fft_notch( dat(:,v), head(v).rate, preproc{ind}.(fnames{k}).freqs );
                end;
                
            case {'fftcomb'},               % FFT comb
                freqs1 = preproc{ind}.(fnames{k}).freqs;
                for v = (1:length(chans))
                    freqs = (mean(freqs1):mean(freqs1):head(v).rate/2)'-(freqs1(2)-freqs1(1))/2;
                    freqs = [ freqs freqs+(freqs1(:,2)-freqs1(:,1)) ];
                    freqs(freqs>head(v).rate/2) = head(v).rate/2;
                    dat(:,v) = fft_notch( dat(:,v), head(v).rate, freqs );
                end;
                
            case {'pca'}                    % PCA
                if ( preproc{ind}.(fnames{k}) )
                    [~,dat] = pca( dat );
                end;
                
            case {'fastica'}                % FastICA
                if ( preproc{ind}.(fnames{k}) )
                    if exist('fastica','file')
                        dat = fastica( dat.' ).';
                    else
                        error(' FastICA not found on the Matlab path!');
                    end;
                end;
                
        end;
    end;
end;
