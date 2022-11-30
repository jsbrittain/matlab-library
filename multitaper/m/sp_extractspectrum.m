function spect = sp_extractspectrum( time, sp, baseline_t )

% Parameters
tlim = dsearchn( time', baseline_t' );
spect = zeros(size(sp,1),size(sp,3));

% Extract spectrum
for n = (1:size(sp,3))
    spect(:,n) = mean( sp(:,(tlim(1):tlim(2)),n), 2);
end;
