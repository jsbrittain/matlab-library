function [E,V]=dpss(N, W, method, int, size)

% Syntax: [E,V]=dpss(N, W); [E,V]=dpss(N, W, method, int, size)
% Dpss generates discrete prolate spheroidal 
% sequences of length N for window width W.
% W should be 2, 5/2, 3, 7/2, 4, etc.
%
% Returns: E: matrix of the first 2W eigenvectors
%          V: vector of the first 2W eigenvalues.
% 
% If only N and W are given, dpss first looks for a dpss
% file on disk. If this is not found, the dpss are generated
% by spline interpolation from dpss of length 256. If these
% are not found, the dpss are calculated directly.
%
% These defaults can be overridden by the other arguments:
%
% 'Method' should be either 'calc' or 'terp'.
%     'calc': Compute the dpss using dpsscalc.m;
%     'terp': Interpolate new dpss from existing
%             ones using dpssint.m.
% Interpolation is significantly faster than calculation.
%
% If 'terp' is specified, the interpolation method 'int' 
% can be specified as 'linear' or 'spline' (the default).
% Linear is faster, but splining is slightly better.
%
% The size of the dpss to be used in the interpolation can
% be specified in 'size' - these dpss must be stored on disk.
% The default is 256. 
%
% For further details, see dpsscalc.m and dpssint.m.
%
% Written by Eric Breitenberger, version date 10/3/95.
% Please send comments and suggestions to eric@gi.alaska.edu
%

if nargin==2
  if rem(W,1)==0
    filename=['s' int2str(N) '_' int2str(W) '.mat']; 
  else
    filename=['s' int2str(N) '_' int2str(10*W) '.mat']; 
  end

  if exist(filename)==2
    msg=['Loading the DPSS from disk...']; 
    disp(msg)
    eval(['load ' filename])
  else
    if rem(W,1)==0
      filename=['s' int2str(256) '_' int2str(W) '.mat']; 
    else
      filename=['s' int2str(256) '_' int2str(10*W) '.mat']; 
    end
    if exist(filename)==2
      disp('Interpolating the DPSS ...')
      [E,V]=dpssint(N,W);
    else
      disp('Calculating the DPSS ...')
      [E,V]=dpsscalc(N,W);
    end
  end

elseif nargin==3 
  if isstr(method)
    if length(method)==6
      msg=['Computing the DPSS using ' method ' interpolation...']; 
      disp(msg)
      [E,V]=dpssint(N,W,method);
    elseif length(method)==4
      if method=='calc'
        msg=['Computing the DPSS using direct calculation...']; 
        disp(msg)
        [E,V]=dpsscalc(N,W);
      elseif method=='terp'
        disp('Interpolating the DPSS ...')
        [E,V]=dpssint(N,W);
      else
        error('Improper specification of arguments.')
      end
    end
  else
    msg=['Computing DPSS using spline interpolation, size ' int2str(method) '...'];
    disp(msg)
    [E,V]=dpssint(N,W,method);
  end

elseif nargin==4
  if length(method)==6
    msg=['Computing the DPSS using ' method ' interpolation, size ' int2str(int) '...']; 
    disp(msg)
    [E,V]=dpssint(N,W,int,method);
  else
    error('Improper specification of arguments.')
  end

else
  msg=['Computing the DPSS using ' int ' interpolation, size ' int2str(size) '...']; 
  disp(msg)
  [E,V]=dpssint(N,W,size,int);
end

k=2*W;  % Return only first k values
E=E(:,1:k);
V=V(1:k);

