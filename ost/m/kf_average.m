function [spavg,paramsavg] = kf_average( sp, params )

if ((~iscell(sp)) | (~iscell(params)))
    error(' This routine is designed to average and scale OST estimates appropriately.');
end;
if ((length(sp) ~= length(params)) | (length(sp)<1))
    error(' Empty cells or input parameters of differing size.');
end;
if (params{1}.jackknife)
    error(' Haven''t thought about whether this routine is valid for jackknife estimates yet!');
end;
if (params{1}.method~=2)
    warning(' Haven`t checked output on anything other than k2 method yet!');
end;

% average trigger zones
        
M = length(sp);
spm = sp{1}; paramsavg = params{1};
for n=(2:M)
    spm = cat(4,spm,sp{n});
    paramsavg.Psp  = paramsavg.Psp + params{n}.Psp;     % Add variances for average
end;
paramsavg.Psp = paramsavg.Psp / M^2;                    % Scale summed variances

spavg(:,1,:) = exp(mean(log(spm(:,1,:,:)),4));
spavg(:,2,:) = exp(mean(log(spm(:,2,:,:)),4));
spavg(:,3,:) = real(exp(mean(log(real(spm(:,3,:,:))),4))) + i*real(exp(mean(log(imag(spm(:,3,:,:))),4)));
spavg(:,4,:) = abs(spavg(:,3,:)).^2./(spavg(:,1,:).*spavg(:,2,:));
spavg(:,5,:) = angle(spavg(:,3,:));

paramsavg.L = paramsavg.L * M;
paramsavg.Lprime = paramsavg.Lprime * M;
