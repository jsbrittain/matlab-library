function varargin = trans( varargin )

for n = (1:nargin)
    varargin{n} = varargin{n}';
end;
