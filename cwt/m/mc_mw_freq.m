function mc=mc_mw_freq(trials,secs,dt,df,frange,opt,mother,wlopt);
%function mc=mc_mw_freq(trials,secs,dt,df,frange,opt,mother,wlopt);
%
% Function to perform monte-carlo simulations using scale oriented wavelet
% transforms. This function performs monte-carlo simulations on simulated
% gaussian white noise signals and returns the 95% confidence limit of
% single trial wavelet coherence between two such such signals obtained
% through Heisenberg box smoothing on the Time-Frequency plane.
%
% EVALUATES COHERENCE ONE SCALE AT A TIME FOR MEMORY CONSERVATION
% REQUIRES FRANGE TO BE SPECIFIED AS [FMIN FMAX], NO OTHER FORMATS ARE
% ACCEPTED AT THIS TIME.
%
%function mc=mc_mw_freq(trials,secs,dt,df,frange,opt,mother,wlopt);

% Check frange structure - SEE ABOVE
if (length(frange)~=2)
    error(' Multi-wavelet monte-carlo simulation failed: FRANGE must be [fmin fmax] at this time');
end;
f_min=frange(1);
f_max=frange(2);

% Data parameters
rate=1/dt;
N=secs*rate;

maxk=200;                   % Arbitrary upper limit for max K
zeta=0.95;                  % Energy concentration to determine max K
be=wlopt{1}; ga=wlopt{2}; A=wlopt{3}; r=(2*be+1)/ga; % wlopt={be,ga,A}
C=ga*A*gamma(r)^2/(gamma(r+1-1/ga)*gamma(r+1/ga))+1;
E=betainc((C-1)/(C+1),[0:maxk-1]+1,r-1).^2;
maxk=max(find(E>=zeta));
dk=E(1:maxk)';              % Energy [0,1] used as weight
weights=dk./sum(dk);

% Calculate scales
wlparam=wlp_morse(wlopt{:},dk);%,wlopt{end});     % Pass full multi-wavelet params {beta,gamma,A,dk}
s_max=1./(wlparam.fourier_factor*f_min);
s_min=1./(wlparam.fourier_factor*f_max);
% Convert min frequency to max scale
s0=s_min;
dj=df;%/wlparam.fourier_factor;
J=ceil((1/dj)*log2(s_max/s0));
wlparam.scale=s0.*2.^([0:J]*dj);
sizeW=[length(wlparam.scale) N]
% Determine numerical COI
wlparam.dt=dt;
wlparam.dj=dj;
wlparam.srange=[s0 J];
wlparam.opt=opt;
wlparam.wlopt=wlopt;
wlparam.multiwavelet=logical(1);
wlparam=numericalCOI(wlparam,N);

% Determine COI region
nbins=1000;
outsidecoi=logical(zeros(sizeW));
for ind=1:length(wlparam.scale)
    outsidecoi(ind,:)=(wlparam.scale(ind)<=wlparam.coi);
end;
maxscale=1;
for ind=1:length(wlparam.scale)
    if any(outsidecoi(ind,:)>0)
        maxscale=ind;
    else
        sig95(ind)=NaN;
    end;
end;
wlc=zeros(maxscale,nbins);
scales=wlparam.scale;
scale_at_a_time=logical(1);

if (scale_at_a_time)    
    wlf_mother=['wlf_' mother];
    halflength=floor(N/2);
	omega=[(2*pi*(0:halflength)/(N*dt)) (-(2*pi*((halflength+1):(N-1)))/(N*dt))]';
    
    wlf=zeros(N,maxk,maxscale);
    for k=1:maxk
        wlopt1=wlopt; wlopt1{end}=k-1; wlopt1{end+1}=dk;
        for ind=1:maxscale
            wlf(:,k,ind)=conj(feval(wlf_mother,omega,scales(ind),dt,wlopt1{:}));
        end;
    end;
end;

% Perform Monte-Carlo analysis
wb=waitbar(0,'Computing monte-carlo simulations...','Name','Multi-Wavelet coherence');
for trial=1:trials
    disp(['Trial ' int2str(trial) ' of ' int2str(trials)]);
    
    % Generate data
    dat1=randn(N,1);
    dat2=randn(N,1);
    
    if (scale_at_a_time)    % Transform one scale at a time - Better memory requirements
        fx1=fft(dat1); fx2=fft(dat2);
        % Traverse scales one at a time (memory saving)
        for ind=1:maxscale
            % Update display
            disp(['Trial ' int2str(trial) ', Scale ' int2str(ind) ' of ' int2str(maxscale)]);
            S11=zeros(1,N); S22=zeros(1,N); S12=zeros(1,N);
            for k=1:maxk
                % Compute normalised wavelet transform at scale s
                W1=ifft(fx1.*wlf(:,k,ind))';
                W2=ifft(fx2.*wlf(:,k,ind))';
                
                % Recursively form multi-wavelet spectra
                S11=S11+weights(k)*(abs(W1).^2);
                S22=S22+weights(k)*(abs(W2).^2);
                S12=S12+weights(k)*(W1.*conj(W2));
            end;
            % Generate coherence
            warning off; WCo=abs(S12).^2./(S11.*S22); warning on;
            clear('S11','S22','S12');
            % Confidence limits (See Grinsted routines)
            cd=WCo(find(outsidecoi(ind,:)));
            cd=max(min(cd,1),0);                    % Ensure normalised [0,1]
            cd=floor(cd*(nbins-1))+1;
            for ind2=1:length(cd)
                wlc(ind,cd(ind2))=wlc(ind,cd(ind2))+1;
            end;
        end;
    else                    % Use full transform - Must quicker if possible
        % Wavelet transform
        [S11,S22,S12,wlparam]=mwTransform(dat1,dat2,dt,dj,scales(1:maxscale),opt,mother,wlopt);       
        % Generate coherence
        warning off; WCo=abs(S12).^2./(S11.*S22); warning on;
        clear('S11','S22','S12');
        % Recurse scales one at a time
        for ind=1:maxscale
            % Confidence limits (See Grinsted routines)
            cd=WCo(ind,find(outsidecoi(ind,:)));
            cd=max(min(cd,1),0);                    % Ensure normalised [0,1]
            cd=floor(cd*(nbins-1))+1;
            for ind2=1:length(cd)
                wlc(ind,cd(ind2))=wlc(ind,cd(ind2))+1;
            end;
        end;
    end;
    
    % Update progress bar
    waitbar(trial/trials);
end;
wlparam.scale=scales;
close(wb);

% Determine 95% confidence limit (See Grinsted routines)
for ind=1:maxscale
    % Bias
    rsqy=((1:nbins)-.5)/nbins;
    bias(ind)=sum(wlc(ind,:).*rsqy)/sum(wlc(ind,:));
    % Variance
    variance(ind)=sum(wlc(ind,:).*(rsqy.^2))/sum(wlc(ind,:))-bias(ind)^2;
    % Confidence limit
    rsqy=((1:nbins)-.5)/nbins;
    ptile=wlc(ind,:);
    idx=find(ptile~=0);
    ptile=ptile(idx);
    rsqy=rsqy(idx);
    ptile=cumsum(ptile);
    ptile=(ptile-.5)/ptile(end);
    if (length(ptile)==1)           % All values fall in same bin (can occur for low freqs and trials==1)
        R95(ind)=NaN;
    else
        R95(ind)=interp1(ptile,rsqy,.95);
    end;
end

% Return wlparam with monte-carlo results
mc.wlparam=wlparam;
mc.wave_opt=opt;
mc.trials=trials;
mc.scale=wlparam.scale(1:maxscale);
mc.freqs=1./(wlparam.fourier_factor*mc.scale);
mc.wlc=wlc;
mc.nbins=nbins;
mc.maxscale=maxscale;
mc.bias=bias;
mc.variance=variance;
mc.R95=R95;
