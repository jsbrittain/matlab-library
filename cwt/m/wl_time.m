function time = wl_time(time,S)
%function time = wl_time(time,S)

time = time(1) + ((0:(size(S,2)-1))/(size(S,2)-1))*(time(end)-time(1));
