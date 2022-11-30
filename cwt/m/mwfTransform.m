function [Sx,out2,Sxy,wlparam]=mwfTransform(varargin)
%function [Sx[,Sy,Sxy],wlparam]=mwfTransform(x[,y],dt,df[,frange],opt,mother,wlopt);
%
% Multiwavelet (MW) Transform
%
% Input parameters
%   x           Time-series ch.1
%   y           (opt) Time-series ch.2
%   dt          Time resolution
%   df          Frequency resolution 1/(scales/octave)
%   frange      (opt) Frequency range [fmin fmax] ([]=default)
%   opt         Analysis options (as wlfTransform)
%                   Extras
%                       e   Return eigenspectra
%   mother      Mother wavelet
%   wlopt       Wavelet options
%
% Weighting schemes (internal variable)
%   0 - Uniform
%   1 - Eigenspectra (default)
%
% Supported wavelet
%   'morse'     {beta, gamma, A(area)}
%               Negative `A' forces that number of taper, with uniform weighting
%
%function [Sx[,Sy,Sxy],wlparam]=mwfTransform(x[,y],dt,df[,frange],opt,mother,wlopt);

% Updated 24/02/2006

% Routine parameters
weighting_scheme=1;         % 0-Unbiased; 1-Eigenspectra
maxk=200;                   % Arbitrary upper limit for max K
eigenspectra = false;
zeta = 0.95;

% Check number of input time series
if ~(length(varargin{1})>1)                     % Length of 'x'
    error(' Time series too short (Must contain multiple elements).');
end;
if (length(varargin{2})==length(varargin{1}))   % Two time series (len 'y' == len 'x')
    bivariate=true;
else                                            % One time series
    bivariate=false;
end;

% Determine input arguments
if (bivariate)
    varargin1={varargin{[1 3:end]}};
    varargin2={varargin{[2:end]}};
else
    varargin1=varargin;
end;

% Calculate energy concentrations
wlopt=varargin{end};
mother=varargin{end-1};
switch (mother)
    case 'morse'
        be=wlopt{1}; ga=wlopt{2}; A=wlopt{3}; r=(2*be+1)/ga;        % wlopt={be,ga,A}
        if (A>0)
            C=ga*(2*A)*gamma(r)^2/(gamma(r+1-1/ga)*gamma(r+1/ga))+1;    % NB: A=|D|/2
            E=betainc((C-1)/(C+1),[0:maxk-1]+1,r-1).^2;
            maxk=max(find(E>=zeta));
        else
            maxk = -A;
            weighting_scheme = 0;
        end;
    otherwise
        error([' Unknown multiwavelet ''' mother '']);
end;

% Determine eigenspectra weighting
switch (weighting_scheme)
    case 0              % Uniform
        dk=ones(maxk,1);
    case 1              % Eigenspectra based
        dk=E(1:maxk).';
end;
weights=dk./sum(dk);

% Perform wavelet transform at each wavelet order
for k = (1:maxk)
    % Display progress
    disp(['Multi-wavelet transform ' int2str(k) ' of ' int2str(maxk)]);
    
    % Perform transform for kth-order wavelet
    argin1=varargin1; argin1{end}{end}=k-1; argin1{end}{end+1}=dk;         % wlopt={...,k,dk}
    [W1,wlparam]=wlfTransform(argin1{:});
    if (bivariate)
        argin2=varargin2; argin2{end}{end}=k-1; argin2{end}{end+1}=dk;    % (as above)
        [W2,wlparam]=wlfTransform(argin2{:});
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
    
    % Record eigenspectra
    if (eigenspectra)
        if (k==1)
            wlparam.eigenspectra.Sx = zeros(size(Sx,1),size(Sx,2),maxk);
            if (bivariate)
                wlparam.eigenspectra.Sy = zeros(size(Sx,1),size(Sx,2),maxk);
                wlparam.eigenspectra.Sxy = zeros(size(Sx,1),size(Sx,2),maxk);
            end;
        end;
        wlparam.eigenspectra.Sx(:,:,k) = (abs(W1).^2);
        if (bivariate)
            wlparam.eigenspectra.Sy(:,:,k) = (abs(W2).^2);
            wlparam.eigenspectra.Sxy(:,:,k) = (W1.*conj(W2));
        end;
    end;
end;

% Define parameters structure
wlparam.dk=dk;
wlparam.L=1/sum((dk/sum(dk)).^2);       % Effective No. segments
wlparam.wlopt=varargin{end};
wlparam.multiwavelet=true;

% Define output argments
if (bivariate)
    out2=Sy;
else
    out2=wlparam;
    clear('wlparam');
end;
