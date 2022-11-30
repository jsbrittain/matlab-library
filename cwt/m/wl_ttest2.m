function p = wl_ttest2( spArray1, spArray2, paired, logtransform )
%
% Takes two 3D array spArray of ( frequency x time x reps )
%  and performs two-scample t-tests on the each time-frequency point, with
%  the null hypothesis of equal mean.
%
%

% Determine input parameters
if (~exist('logtransform'))
    logtransform = [];
end;
if (~exist('paired'))
    paired = [];
end;

% Defaults
if (isempty(logtransform))
    logtransform = true;
end;
if (isempty(paired))
    paired = true;
end;

% Variance-normalise spectra
if (logtransform)
    disp('Log-transforming spectral data...');
    spArray1 = log( spArray1 ); spArray1(isinf( spArray1 )) = 1/eps;
    spArray2 = log( spArray2 ); spArray1(isinf( spArray2 )) = 1/eps;
end;

% Pre-allocate memory
p = zeros(size(spArray1,1),size(spArray1,2));
disp('Computing t-test on spectrograms...');

% Perform T-test per (time,frequency) tuple
for x = (1:size(spArray1,2))
    for y = (1:size(spArray1,1))
        if (paired)
            [h,p(y,x)] = ttest(squeeze(spArray1(y,x,:)),squeeze(spArray2(y,x,:)));
        else
            [h,p(y,x)] = ttest2(squeeze(spArray1(y,x,:)),squeeze(spArray2(y,x,:)));
        end;
    end;
end;

% Report progress
disp('done.');
