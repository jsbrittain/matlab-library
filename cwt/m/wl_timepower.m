function timecourse = wl_timepower(S,params,frange)
%function timecourse = wl_timepower(S,params,frange)

fsearch = dsearchn(params.freqs',frange');
timecourse = mean(S(fsearch(1):fsearch(2),:),1);
