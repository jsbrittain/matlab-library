function psp_ch1(f,cl,freq,label)
% function to plot coherence in current figure/subplot window
%  psp_chi(f,cl,freq,label)
%
% f,cl     Output from pooled spectral analysis routine.
% freq     Frequency limit for plotting (Hz).
% label    Optional title instead of cl.what.

freq_pts=round(freq/cl.df);
%Check freq range
[x,y]=size(f);
if (freq_pts>x)
	error('Requested frequency range too large.');
end

f_max=freq_pts*cl.df;
chi_cl=[cl.df,cl.chi_c95;f_max,cl.chi_c95];
plot(f(1:freq_pts,1),f(1:freq_pts,8),'k-',chi_cl(:,1),chi_cl(:,2),'k--');
axis([0,freq,-Inf,Inf]);
if (nargin>4)
  title(label);
else
 title(['\chi^2 test: ',cl.what]);
end
