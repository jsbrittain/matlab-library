function [dat,rate]=smr_isolatedat(smrdat,head,channels)

% Check channels
latest_start=0;
dt=channels{1}.sampleinterval;
for ind=(1:length(channels))
    
    if (channels{ind}.sampleinterval~=dt)
        error(' Sampling intervals inconsistent (channel ' num2str(ind) ').');
    end;
    if (~isfield(channels{ind},'start'))
        error(' Channel ' num2str(ind) ' is not a data series.');
    end;
    
    latest_start = max([latest_start head{ind}.start]);
    earliest_stop = min([earliest_stop head{ind}.stop]);
    
end;

% Reform data
for ind=(1:length(channels))
    out(:,ind) = dat{ind}(tt);
end;
