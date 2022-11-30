function [Sx,out2,Sxy,wlparam]=mwTransform(varargin)
%function [Sx[,Sy,Sxy],wlparam]=mwTransform(x[,y],dt,dj[,srange],opt,mother,wlopt);
%
% Multiwavelet (MW) Transform
%
% Input parameters
%   x           Time-series ch.1
%   y           (opt) Time-series ch.2
%   dt          Time resolution
%   dj          Scale resolution
%   srange      (opt) Scale range [s0 J] (gives s0*2^([0:J]*dj))
%   opt         Analysis options (as wlTransform)
%   mother      Mother wavelet
%   wlopt       Wavelet options
%
% Weighting schemes (internal variable)
%   0 - Uniform
%   1 - Eigenspectra (default)
%
% Supported wavelet
%   'morse'     {beta, gamma, A(area)}
%
% Weighting schemes
%   0 - Unbiased
%   1 - Eigenspectra (default)
%
%function [Sx[,Sy,Sxy],wlparam]=mwTransform(x[,y],dt,dj[,srange],opt,mother,wlopt);

% Updated 22/11/2005

% Routine parameters
weighting_scheme=1;         % 0-Unbiased; 1-Eigenspectra
maxk=200;                   % Arbitrary upper limit for max K
zeta=0.95;                  % Energy concentration to determine max K

% Check number of input time series
if ~(length(varargin{1})>1)                     % Length of 'x'
    error(' Time series too short (Must contain multiple elements).');
end;
if (length(varargin{2})==length(varargin{1}))   % Two time series (len 'y' == len 'x')
    bivariate=logical(1);
else                                            % One time series
    bivariate=logical(0);
end;

% Determine input arguments
if (bivariate)
    varargin1={varargin{[1 3:end]}};
    varargin2={varargin{[2:end]}};
else
    varargin1=varargin;
end;
wlopt=varargin{end};
be=wlopt{1}; ga=wlopt{2}; A=wlopt{3}; r=(2*be+1)/ga; % wlopt={be,ga,A}
C=ga*(A)*gamma(r)^2/(gamma(r+1-1/ga)*gamma(r+1/ga))+1; % NB: A=|D|/2
E=betainc((C-1)/(C+1),[0:maxk-1]+1,r-1).^2;
maxk=max(find(E>=zeta));

% Determine eigenspectra weighting
switch (weighting_scheme)
    case 0              % Unbiased
        dk=ones(maxk,1);
    case 1              % Eigenspectra based
        dk=E(1:maxk)';  % Energy [0,1] used as weight
end;
weights=dk./sum(dk);

% Perform wavelet transform at each wavelet order
for k=1:maxk
    % Perform transform for kth-order wavelet
    argin1=varargin1; argin1{end}{end}=k-1; argin1{end}{end+1}=dk;         % wlopt={...,k,dk}
    [W1,wlparam]=wlTransform(argin1{:});
    if (bivariate)
        argin2=varargin2; argin2{end}{end}=k-1; argin2{end}{end+1}=dk;    % (as above)
        [W2,wlparam]=wlTransform(argin2{:});
    end;
    
    % Initialise recusive variables
    if (k==1)
        Sx=zeros(size(W1));
        if (bivariate)
            Sy=zeros(size(W1));
            Sxy=zeros(size(W1));
        end;
    end;
    
    % Recursively form multi-wavelet spectra
    Sx=Sx+weights(k)*(abs(W1).^2);
    if (bivariate)
        Sy =Sy +weights(k)*(abs(W2).^2);
        Sxy=Sxy+weights(k)*(W1.*conj(W2));
    end;
end;

% Define parameters structure
wlparam.dk=dk;
wlparam.L=1/sum((dk/sum(dk)).^2);       % Effective No. segments
wlparam.wlopt=wlopt;

% Define output argments
if (bivariate)
    out2=Sy;
else
    out2=wlparam;
    clear('wlparam');
end;
