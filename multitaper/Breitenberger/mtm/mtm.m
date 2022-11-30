 function [S,c]=mtm(x,W, method, conf);
% Syntax: [S,c]=mtm(x,W); [S,C]=mtm(x,W,method,conf);
% Mtm is a shell function, which implements the Thomson
% multiple-taper method of spectral estimation for the
% time series x, for window width W. W should be specified
% as e. g. 2, 5/2, 3, 7/2, etc. W=4 is often a good choice
% for the bandwidth - this is the default.
%
% There are two options for the estimation method ("method"):
%   'lin'  : simple estimate with linear weights.
%   'dap' : adaptively-weighted estimate. (default)
%
% The confidence level for the confidence limits can be
% specified in 'conf'. The default is .95 (95%).
%
% Output: The spectral estimate is returned in S,
%         and confidence limits are returned in c.
%
% For more details, see dpss.m, dpssint.m, dpsscalc.m,
% mtmlin.m, and mtmadap.m
%
% Written by Eric Breitenberger, version date 10/1/95.
% Please send comments and suggestions to eric@gi.alaska.edu
%

if nargin==1, W=4; method='dap'; conf=.95;
elseif nargin==2
  if isstr(W), method=W; W=4; conf=.95;
  elseif W<1, method='dap'; conf=W; W=4;
  else method='dap', conf=.95;
  end
elseif nargin==3, 
  if ~isstr(method)
    conf=method; method='dap';
  else, conf=.95;
  end
end

% Get the dpss, one way or another:
N=length(x);
[E,V]=dpss(N,W);

if method=='lin',
  [S,c]=mtmlin(x,E,V);
elseif method=='dap',
  [S,c]=mtmadap(x,E,V);
else
  error('Improper specification of method.')
end

