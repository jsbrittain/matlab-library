
addpath('~/matlab/library/hmm/m');
addpath('~/matlab/library/thirdparty/arfit/arfit');

%% Generate data

switch (2)
    case 0,
        
        m=1; M=2; p=3;
        N=200;
        
        y0={}; A={}; k=1;
        C = rand(1)*rand(m).*eye(m) + 1e-4*rand(1)*rand(m);
        lastwarn('');
        while (k<=M)
            A{k} = rand(m,m*p);
            y0{k} = arsim(zeros(1,m),A{k},C,N).';
            if (isempty(lastwarn))      % Ensure AR model stable (warning output from arsim)
                k=k+1;
            else
                lastwarn('');
            end;
        end;
        y = [y0{:}];
        y = y-repmat(mean(y,2),1,size(y,2));
    case 1,
        
        rate=100;
        N=3*rate;
        p=8;
        
        s2static=0.1^2;
        et=sqrt(s2static)*randn(1,N);
        segN = N/3;
        
        f = [10 20 30];
        y = [sin(2*pi*f(1)*[0:(segN-1)]/rate) sin(2*pi*f(2)*[0:(segN-1)]/rate) sin(2*pi*f(3)*[0:(segN-1)]/rate)] + et;
        state = [1*ones(1,segN) 2*ones(1,segN) 3*ones(1,segN)];
        
    case 2,
        
        % Offset double sinusoids
        
        rate=100;
        N=3*rate;
        m=3; p=8; M=1;
        
        rot = randi(4,1,m)+3;
        rot(1)=0;
        
        s2static=0.1^2; y=[];
        for k=(1:m)
            et=sqrt(s2static)*randn(1,N);
            segN = N/3;

            f = [10 20 30];
            y(k,:) = [sin(2*pi*f(1)*[0:(segN-1)]/rate) sin(2*pi*f(2)*[0:(segN-1)]/rate) sin(2*pi*f(3)*[0:(segN-1)]/rate)] + et;
            y(k,:) = [y(k,(1+rot(k)):end) y(k,1:rot(k))];
        end;
end;

%% VAR analysis
s2static=[]; alpha=0.4;
[q,s2,evidence,ct,rhoavg] = var_kalman(y,p,s2static,alpha);

% Pad states if AR model (AR->VAR)
m2p=m^2*p;
if (length(size(q))<4)
    q = reshape(q,1,1,size(q,1),size(q,2));
    ct = reshape(ct,1,1,size(ct,2));
    m2p=m;
end;

% Plot VAR results
figure;
subplot(p+3,1,1); plot(log(evidence));
subplot(p+3,1,2); plot(y');
subplot(p+3,1,3); hold('on');
strcol='bgrcmyk';
    for n=(1:m2p)
        plot(squeeze(ct(n,n,:)),strcol(mod(n-1,length(strcol))+1));
    end;
for k=(1:p)
    subplot(p+3,1,k+3); hold('on');
    for n1=(1:m)
        for n2=(1:m)
            plot(squeeze(q(n1,n2,k,:)),strcol(mod((n1-1)*m+n2-1,length(strcol))+1));
        end;
    end;
end;
for k=(1:(p+3))
    subplot(p+3,1,k); hold('on');
    for n=(1:(M-1))
        plot(N*n*[1 1],ylim,'k');
    end;
end;
