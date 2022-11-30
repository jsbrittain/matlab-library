function [epochs,erptime]=epochdata_disjoint(dat,rate,duration,baseline,opt_str)
%function [epochs,erptime]=epochdata_disjoint(dat,rate,duration,baseline,opt_str)
%
%
%
%function [epochs,erptime]=epochdata_disjoint(dat,rate,duration,baseline,opt_str)

offset = 0;
trig = (1:duration*rate/1000:(size(dat,1)-duration*rate/1000));
[epochs,erptime]=epochdata(dat,trig,rate,offset,duration,baseline,opt_str);
