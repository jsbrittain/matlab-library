function y=betainc_cutoff(x,a,b,cutoff);
%function y=betainc_cutoff(x,a,b,cutoff);
%
% Returns the absolute difference between the beta function with parameters
% (x,a,b) and the cutoff which lies in the range [0,1].  Used in conjunction
% with a function minimisation routine to determine the inverse incomplete
% beta function.
%
%function y=betainc_cutoff(x,a,b,cutoff);

y=abs(betainc(x,a,b)-cutoff);
