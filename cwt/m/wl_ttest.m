function p = wl_ttest( spArray, logtransform, test )
%
% Take a 3D array spArray of ( frequency x time x reps )
%  and perform t-tests on the each time-frequency point, with the
%  null hypothesis of a zero mean.
%
%

if (~exist('logtransform'))
    logtransform = [];
end;
if (~exist('test'))
    test = [];
end;

if (isempty(logtransform))
    logtransform = true;
end;
if (isempty(test))
    test = 'ttest';
end;

if (logtransform)
    disp('Log-transforming spectral data...');
    spArray = log( spArray );
end;

p = zeros(size(spArray,1),size(spArray,2));
disp('Computing t-test on spectrograms...');

for x = (1:size(spArray,2))
    for y = (1:size(spArray,1))
        switch ( test )
            case 'ttest',     % T-test
                [~,p(y,x)] = ttest(squeeze(spArray(y,x,:)));
            case 'signrank',     % Sign rank (Wilcoxon)
                [p(y,x),~] = signrank(squeeze(spArray(y,x,:)));
            otherwise
                error([' Unknown test - ' test ]);
        end;
    end;
end;

disp('done.');
