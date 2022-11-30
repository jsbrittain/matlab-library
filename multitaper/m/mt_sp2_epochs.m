function [sp11,sp22,sp12,params]=mt_sp2_epochs(dat1,dat2,rate,opt_str)
%function [sp11,sp22,sp12,params]=mt_sp2_epochs(dat1,dat2,rate,opt_str)
%
% Multi-taper analysis type 2
%
% See mt_sp2
%
%function [sp11,sp22,sp12,params]=mt_sp2_epochs(dat1,dat2,trig,rate,opt_str)

[sp11,sp22,sp12,params]=mt_sp0_epochs(dat1,dat2,rate,opt_str);
