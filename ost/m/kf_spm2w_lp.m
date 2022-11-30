function [sp,params,spk] = kf_spm2w_lp(dat1,dat2,trig,offset,duration,rate,seg_pwr,opt_str)
%function [sp,params,spk] = kf_spm2w_lp(dat1,dat2,trig,offset,duration,rate,seg_pwr,opt_str)
%
%
%
% Perform log-filtering but apply weights directly in the linear plane (lp)
%
%function [sp,params,spk] = kf_spm2w_lp(dat1,dat2,trig,offset,duration,rate,seg_pwr,opt_str)

% Default parameters
method = 0;                         % Smoothing method

% Parse opt_str for Kalman-Fourier related options
options=deblank(opt_str);
opt_str='';
while (any(options))                % Parse individual options from string.
	[opt,options]=strtok(options);
	optarg=opt(2:length(opt));      % Determine option argument.
	switch (opt(1))
        case 'k'                    % 
            method=str2num(optarg);
            if (isempty(method))
                error([' Error in parameter ' opt]);
            end;
    	otherwise                   % Pass on to analysis
            opt_str=[opt_str ' ' opt];
    end;
end;

% Retrieve unprocessed segments
[spk,params]=kf_spm2w(dat1,dat2,trig,offset,duration,rate,seg_pwr,[opt_str ' k' num2str(method)]);
[sp0,params0]=kf_spm2w(dat1,dat2,trig,offset,duration,rate,seg_pwr,[opt_str ' k0']);

% Weight by Kalman filter (offline)
sp = zeros(size(spk,1),3,size(spk,3));
M = size(params.w,1); sp0 = sp0(:,(1:3),:);

progress = 0;
disp(' Post-processing weights ...');
for ind=(1:M)
    
    % Display progress
    if (floor(ind/M*100)>progress)
        progress=floor(ind/M*100);
        disp(['  weighting segments: ' int2str(progress) '%']);
    end;
    
    % Form linear plane weighting
    sp(:,:,ind) = sum( sp0 .* reshape( params.w(ind,:,:), size(sp0,1), 3, size(sp0,3) ), 3);
    
end;
disp(' done.');

% Coherence
warning off
sp(:,4,:) = abs(sp(:,3,:)).^2./(sp(:,1,:).*sp(:,2,:));
warning on

% Phase
sp(:,5,:) = angle(sp(:,3,:));
