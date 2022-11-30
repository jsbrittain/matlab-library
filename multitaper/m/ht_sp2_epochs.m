function [sp11,sp22,sp12,params]=ht_sp2_epochs(dat1,dat2,rate,opt_str)
%function [sp11,sp22,sp12,params]=ht_sp2_epochs(dat1,dat2,rate,opt_str)
%
% Multitaper time-frequency analysis utilising Hermite functions
% Spectral analysis over trials
%
% Data input: Epoched
%
% Implementation
%   Eigenvalue weighting
%
% Input parameters
%       dat1        Epoched time series ch.1
%       dat2        Epoched time series ch.2
%       rate        Sampling rate
%       opt_str     Options string
%                       r<n>        rectify channels (n={0:ch.1,1:ch.2,2:ch1,2})
%                       f<n>        maximum analysis frequency
%                       A<n>        localisation area
%                       t<msecs>    time skip for plotting (default: 50 msecs)
%                       s<n>        time scaling for elliptical discs (default: 0.010)
%                       n           power normalise each increment
%       Special options:
%                       e           return eigenspectra matrix
%                                       analysis centred on time mid-point
%                                       index 't' corresponds to eigenspectra
%
%function [sp11,sp22,sp12,params]=ht_sp2_epochs(dat1,dat2,rate,opt_str)

% Squeeze data matrices (converts EEGLab structure to 2D matrix if required)
dat1=squeeze(dat1);
dat2=squeeze(dat2);

% Default parameters
return_eigenspectra = false;

% Parse options string
options=deblank(opt_str);
opt_str2 = '';
while (any(options))                % Parse individual options from string.
	[opt,options]=strtok(options);
	optarg=opt(2:length(opt));      % Determine option argument.
	switch (opt(1))
	    case 'e'                    % Return eigenspectrograms
        	return_eigenspectra = true;
        otherwise
            opt_str2 = [ opt_str2 ' ' opt ];
    end;
end;

% Perform hermite TFR analysis over epochs
trigcount=size(dat1,2);
for ind=(1:trigcount)
    % Display progress
    disp(['Trial ' int2str(ind) ' of ' int2str(trigcount)]);
    % Perform single epoch analysis
    [p11,p22,p12,params]=ht_sp(dat1(:,ind),dat2(:,ind),rate,opt_str2);
    % Form recursive estimate of the spectra
    if ( ind == 1 )
        if ( return_eigenspectra )
            sp11 = zeros([size(p11) trigcount],'single');
            sp22 = sp11; sp12 = sp11;
            sp11(:,:,1)=p11;
            sp22(:,:,1)=p22;
            sp12(:,:,1)=p12;
        else
            sp11=p11/trigcount;
            sp22=p22/trigcount;
            sp12=p12/trigcount;
        end;
    else
        if ( return_eigenspectra )
            sp11(:,:,ind) = p11;
            sp22(:,:,ind) = p22;
            sp12(:,:,ind) = p12;
        else
            % Average as you go
            sp11=sp11+p11/trigcount;
            sp22=sp22+p22/trigcount;
            sp12=sp12+p12/trigcount;
        end;
    end;
end;

% Update parameters structure
params.trigcount=trigcount*params.L;        % For conf limits
