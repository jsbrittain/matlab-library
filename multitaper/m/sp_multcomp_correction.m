function P = sp_multcomp_correction( sp, method, opt_str )
%function P = sp_multcomp_correction( p, method, opt_str )
%
% Multiple comparions correction, taking a permuted probability
% matrix and returning significance zones that survive correction.
%
% Input parameters
%       p           Probability matrix containing permutations.
%                   Must be of form ( x, y, iterations + 1 ), x or y may be unitary.
%                   The first element (iteration=1) is the original (unpermuted) data
%       method      (opt) Inference method
%                                 0 - Null
%                                 1 - Single threshold test (never seen this work properly - is it coded correctly?)
%                                 2 - Maximal suprathreshold cluster size
%                       (default) 3 - Exceedance mass
%       opt_str     (opt) Options string
%               
% To be used with:
%   wl_ttest2_permute.m
%
% Reference
%   Nichols, T. E. and Holmes, A. P. (2001) Nonparametric tests for
%   functional neuroimaging: a primer with examples. Human Brain Mapping
%   15:1-25.
%
%function P = sp_multcomp_correction( p, opt_str )

% Check input parameters
if (~exist('opt_str','var'))
    opt_str = '';
end;
if (~exist('method','var'))
    method = [];
end;

% Defaults for input parameters
if (isempty(method))
    method = 3;
end;

% Option defaults
alpha = 0.05;           % Significance threshold
connected = 4;          % 4-connected or 8-connected zoning method

% Parse options string
options=deblank(opt_str);
while (any(options))                % Parse individual options from string.
	[opt,options]=strtok(options);
	optarg=opt(2:length(opt));      % Determine option argument.
	switch (opt(1))
	    case 'c'                    % Zonal connection method
        	n=str2double(optarg);
            if ((n~=4) && (n~=8))
                error([' Invalid argument -- ' opt ', zonal connectivity method must be 4 or 8 (see bwlabel for further details)']);
            end;
            connected = n;
        otherwise                   % Unknown option
            error(['Unknown option -- ' opt]);
    end;
end;

% Determine data parameters
iters = size(sp,3)-1;

% Select inference method
switch ( method )
    
    case 0,     % Uncorrected
        
        P = nan(size(sp,1),size(sp,2));
        for n1 = (1:size(sp,1))
            for n2 = (1:size(sp,2))
                P(n1,n2) = sum( sp(n1,n2,:)>=sp(n1,n2,1) )/size(sp,3);
            end;
        end;

    case 1,     % Single threshold test
        
        tmax = zeros(1,iters);
        for iter = (1:iters+1)
            tmax(iter) = max(max(sp(:,:,iter)));
        end;
        stmax = sort( tmax, 'descend' ); stmaxP = stmax( floor(alpha*(iters+1))+1 );
        P = sp(:,:,1) > stmaxP;

    case 2,     % Maximal suprathreshold cluster size
        
        maxclustersize = zeros( 1, iters );
        for iter = (1:iters+1)
            % Isolate iteration
            spnow = sp(:,:,iter);
            % count blob size
            blobs = bwlabel( spnow < alpha, connected );
            clustersize = zeros(1,max(blobs(:)));
            for k = (1:max(blobs(:)))
                clustersize(k) = sum(blobs(:)==k);
            end;
            if (~isempty( clustersize ))
                maxclustersize(iter) =  max( clustersize );
            else
                maxclustersize(iter) = 0;
            end;
        end;
        sortmaxclustersize = sort(maxclustersize,'descend');
        critical_th = sortmaxclustersize( floor(alpha*(iters+1))+1 );
        % construct P threshold for plotting
        spnow = sp(:,:,1);
        blobs = bwlabel( spnow < alpha, connected );
        for k = (1:max(blobs(:)))
            clustersize = sum(blobs(:)==k);
            if (clustersize<critical_th)
                blobs( blobs==k ) = 0;
            end;
        end;
        P = (blobs>0);
        
    case 3,     % Exceedance mass
        
        maxclustermass = zeros( 1, iters );
        for iter = (1:iters+1)
            % Isolate iteration
            spnow = sp(:,:,iter);
            % compute exceedance mass ( integral of static map over threshold for each sig cluster )
            blobs = bwlabel( spnow < alpha, connected );
            clustermass = zeros(1,max(blobs(:)));
            for k = (1:max(blobs(:)))
                clustermass(k) = sum( alpha - spnow(blobs==k) );
            end;
            if (~isempty( clustermass ))
                maxclustermass(iter) =  max( clustermass );
            else
                maxclustermass(iter) = 0;
            end;
        end;
        sortmaxclustermass = sort(maxclustermass,'descend');
        critical_th = sortmaxclustermass( floor(alpha*(iters+1))+1 );
        % construct P threshold for plotting
        spnow = sp(:,:,1);
        blobs = bwlabel( spnow < alpha, connected );
        for k = (1:max(blobs(:)))
            clustermass = sum( alpha - spnow(blobs==k) );
            if (clustermass<critical_th)
                blobs( blobs==k ) = 0;
            end;
        end;
        P = (blobs>0);
        
end;
