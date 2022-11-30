function times=ht_eeglab_times(sp,EEG)
%
% Return time vector for ht_ and wl*_ functions
%

times=EEG.times(1)+((0:(size(sp,2)-1))/(size(sp,2)-1))*diff(EEG.times([1 end]));
