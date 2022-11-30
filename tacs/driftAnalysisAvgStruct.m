function analysis = driftAnalysisAvgStruct( analysisSeg )
%function analysis = driftAnalysisAvgStruct( analysisSeg )
%
% Accepts cell arrays of "analysis" structures
%
% Averages structures over cell arrays
%
% Analysis structure can be a 1d array of analysis structures (such as may
% be used to index analysis from different channels)
%
%function analysis = driftAnalysisAvgStruct( analysisSeg )

% Sort histogram fieldnames from scalar fields
allfnames = fieldnames(analysisSeg{1}(1));
fnames = allfnames;
for k = (length(fnames):-1:1)
    if ( length(analysisSeg{1}(1).(fnames{k})) ~= analysisSeg{1}(1).bincount )
        fnames = fnames([1:(k-1) (k+1):end]);
    end;
end;
% Place "histcount" last
fnames = {fnames{~strcmp(fnames,'histcount')},fnames{strcmp(fnames,'histcount')}};

% Concatenate histograms
analysisAvg = analysisSeg{1};
for c = (1:length(analysisSeg))
    for ch = (1:length(analysisAvg))
        for g = (1:length(fnames))
            if (strcmp(fnames{g},'histcount'))
                normaliser = 1;
            else
                normaliser = analysisSeg{c}(ch).histcount;
            end;
            if (c==1)
                analysisAvg(ch).(fnames{g}) = normaliser.*analysisSeg{c}(ch).(fnames{g});
            else
                analysisAvg(ch).(fnames{g}) = nansum( [ analysisAvg(ch).(fnames{g}); normaliser.*analysisSeg{c}(ch).(fnames{g}) ], 1 );
            end;
        end;
    end;
end;
% Normalise by histogram count
for ch = (1:length(analysisAvg))
    for g = (1:(length(fnames)-1))
        analysisAvg(ch).(fnames{g}) = analysisAvg(ch).(fnames{g})./analysisAvg(ch).histcount;
    end;
end;
% Replace likelihood with histogram count - can be slightly different (imprecision with large sample counts???)
% Histogram counts look more like concatenated data than concatenated histogram approach anyway
for ch = (1:length(analysisAvg))
    analysisAvg(ch).likeli = analysisAvg(ch).histcount/sum(analysisAvg(ch).histcount)/analysisAvg(ch).bincount;
end;

% Remove all old scalar fields
analysisAvg = rmfield(analysisAvg,setdiff(allfnames,fnames));

% Calculate new summary statistics
for ch = (1:length(analysisAvg))
    analysis(ch) = driftAnalysisSummaryStats( analysisAvg(ch) );
end;
