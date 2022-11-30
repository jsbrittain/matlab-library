function [permstats,channamesPerm,analysis] = analysisDriftPerm_segdata( hdat, hdatSham, ch, channames, permutations, bincount )

if (~exist('bincount','var'))
    bincount = [];
end;
if (isempty(bincount))
    bincount = 20;
end;

% Permutation analysis

if (isfield(ch,'acc')),       ch.accs = ch.acc;         end;
if (isfield(ch,'acc_pca')),   ch.accs_pca = ch.acc_pca; end;
if (isfield(ch,'emg')),       ch.emgs = ch.emg;         end;
if (~isfield(ch,'accs')),     ch.accs = [];             end;
if (~isfield(ch,'accs_pca')), ch.accs_pca = [];         end;
if (~isfield(ch,'emgs')),     ch.emgs = [];             end;
if (~isfield(ch,'sep_trig')), ch.septrig = [];          end;    % Triggered analysis rather than TACS
if (length([ch.septrig ch.tacs])~=1)
    error(['Incorrect number of stimulation channels specified: ' length([ch.septrig ch.tacs]) ' (should be 1).']);
end;
chs = [ ch.accs ch.accs_pca ch.emgs ch.septrig ch.tacs ];       % TACS must be last channel!
channamesPerm = channames(chs);

chcount = length(chs)-1;
hdatBoth = cat(2,hdat(:,:,chs),hdatSham(:,:,chs));
set1len = size(hdat,2);
set2len = size(hdatSham,2);
labels = [ ones(1,set1len) zeros(1,set2len) ];
nulldistr = struct('amp',[],'likeli',[],'meanamp',[]);

% Concatenation method; 1=data, 2=epoched-histograms
switch ( 2 )
    case 1,         %%% Old version --- concatenate data
        for k = (1:(permutations+1))
            disp(sprintf('Permutation %g / %g...',k-1,permutations));
            % Last iteration is not shuffled (putting last keeps analysis in memory too)
            if (k==(permutations+1))
                order = (1:length(labels));
            else
                order = randperm(length(labels));
            end;
            hdatGrp1 = reshape( hdatBoth(:,order(1:set1len),:), size(hdatBoth,1)*set1len, size(hdatBoth,3) );
            hdatGrp2 = reshape( hdatBoth(:,order((set1len+1):end),:), size(hdatBoth,1)*set2len, size(hdatBoth,3) );
            % Analysis
            TACSangGrp1 = angle(hdatGrp1(:,end));
            TACSangGrp2 = angle(hdatGrp2(:,end));
            for ch = (1:chcount)
                analysisGrp1(ch) = driftAnalysis( angle(hdatGrp1(:,ch)), TACSangGrp1, abs(hdatGrp1(:,ch)), bincount );
                analysisGrp2(ch) = driftAnalysis( angle(hdatGrp2(:,ch)), TACSangGrp2, abs(hdatGrp2(:,ch)), bincount );
                % Group PSI difference
                nulldistr.meanamp(ch,k) = nanmean(analysisGrp1(ch).amp) - nanmean(analysisGrp2(ch).amp);                % Mean amplitude difference
                nulldistr.meanfreq(ch,k) = nanmean(analysisGrp1(ch).instfreq) - nanmean(analysisGrp2(ch).instfreq);     % Mean frequency difference
                nulldistr.amp(ch,k) = abs( analysisGrp1(ch).amp_psi ) - abs( analysisGrp2(ch).amp_psi );                % PSI amplitude modulation difference
                nulldistr.likeli(ch,k) = abs( analysisGrp1(ch).likeli_psi ) - abs( analysisGrp2(ch).likeli_psi );       % PSI entrainment modulation difference
            end;
        end;
        fnames = fieldnames(nulldistr);
        for k = (1:length(fnames))
            metric.(fnames{k}) = nulldistr.(fnames{k})(:,end);
            nulldistr.(fnames{k}) = nulldistr.(fnames{k})(:,1:(end-1));

            for ch = (1:size(nulldistr.(fnames{k}),1))
                mu = mean(nulldistr.(fnames{k})(ch,:));
                sd = std(nulldistr.(fnames{k})(ch,:));
                tstat.(fnames{k})(ch) = (metric.(fnames{k})(ch)-mu)/sd;
                % One sided p-value
                p1 = 1-find((metric.(fnames{k})(ch)<=sort(nulldistr.(fnames{k})(ch,:))),1,'first')/permutations;
                if (isempty(p1)), p1 = 1/permutations; end; if (p1>0.5), p1=1-p1; end;
                pval1.(fnames{k})(ch) = p1;
                % Two-sided p-value
                p2 = 1-find((abs(metric.(fnames{k})(ch))<=sort(nulldistr.(fnames{k})(ch,:))),1,'first')/permutations + find((-abs(metric.(fnames{k})(ch))<=sort(nulldistr.(fnames{k})(ch,:))),1,'first')/permutations;
                if (isempty(p2)), p2 = 1/permutations; end;
                pval2.(fnames{k})(ch) = p2;
            end;
        end;
        
    case 2,         %%% New version --- concatenate histograms from each segment
        % Construct drift histograms per epoch
        for seg = (1:size(hdatBoth,2))
            for ch = (1:chcount)
                analysisSeg{seg}(ch) = driftAnalysis( angle(hdatBoth(:,seg,ch)), angle(hdatBoth(:,seg,end)), abs(hdatBoth(:,seg,ch)), bincount );
            end;
        end;
        % Permute epochs, computing stats from concatenated histograms
        fprintf('<'); progress = 0;
        for k = (1:(permutations+1))
            newprogress = floor(100*k/(permutations+1));
            if (( newprogress > progress ) && (progress<100) )
                progress = newprogress;
                if (mod(newprogress,10)==0)
                    fprintf('%g',newprogress);
                else
                    fprintf('.');
                end;
            end;
            % Last iteration is not shuffled (putting last keeps analysis in memory too)
            if (k==(permutations+1))
                order = (1:length(labels));
            else
                order = randperm(length(labels));
            end;
            orderset1 = order(1:set1len);
            orderset2 = order((set1len+1):end);
            % Concatenate histograms from analysis structures
            analysisGrp1 = driftAnalysisAvgStruct( analysisSeg(orderset1) );
            analysisGrp2 = driftAnalysisAvgStruct( analysisSeg(orderset2) );
            for ch = (1:chcount)
                % Group PSI difference
                nulldistr.meanamp(ch,k) = nanmean(analysisGrp1(ch).amp) - nanmean(analysisGrp2(ch).amp);                % Mean amplitude difference
                nulldistr.meanfreq(ch,k) = nanmean(analysisGrp1(ch).instfreq) - nanmean(analysisGrp2(ch).instfreq);     % Mean frequency difference
                nulldistr.amp(ch,k) = abs( analysisGrp1(ch).amp_psi ) - abs( analysisGrp2(ch).amp_psi );                % PSI amplitude modulation difference
                nulldistr.likeli(ch,k) = abs( analysisGrp1(ch).likeli_psi ) - abs( analysisGrp2(ch).likeli_psi );       % PSI entrainment modulation difference
                nulldistr.instfreq(ch,k) = abs( analysisGrp1(ch).instfreq_psi ) - abs( analysisGrp2(ch).instfreq_psi );       % PSI frequency difference
                nulldistr.ampPsi1(ch,k,:) = analysisGrp1(ch).amp_psi;
                nulldistr.ampPsi2(ch,k,:) = analysisGrp2(ch).amp_psi;
                nulldistr.likeliPsi1(ch,k,:) = analysisGrp1(ch).likeli_psi;
                nulldistr.likeliPsi2(ch,k,:) = analysisGrp2(ch).likeli_psi;
            end;
        end;
        fprintf('>\n');
end;

% Calcualte t-stats and p-values
fnames = fieldnames(nulldistr);
for k = (1:length(fnames))
    metric.(fnames{k}) = nulldistr.(fnames{k})(:,end);
    nulldistr.(fnames{k}) = nulldistr.(fnames{k})(:,1:(end-1));
    for ch = (1:size(nulldistr.(fnames{k}),1))
        mu = mean(nulldistr.(fnames{k})(ch,:));
        sd = std(nulldistr.(fnames{k})(ch,:));
        tstat.(fnames{k})(ch) = (metric.(fnames{k})(ch)-mu)/sd;
        try
            [pval1.(fnames{k})(ch),pval2.(fnames{k})(ch)] = analysisDriftPermJointProb( metric.(fnames{k})(ch), nulldistr.(fnames{k})(ch,:) );
        catch
            pval1.(fnames{k})(ch) = NaN;
            pval2.(fnames{k})(ch) = NaN;
        end;
    end;
end;

% Collate permutation statistics
permstats.metric = metric;
permstats.nulldistr = nulldistr;
permstats.tstat = tstat;
permstats.p1 = pval1;
permstats.p2 = pval2;

% Return unpermuted group analyses
analysis = [ analysisGrp1; analysisGrp2 ];
