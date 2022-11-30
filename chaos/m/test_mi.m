
% Time-series
T = 1e6;                            % Works best with large samples!
x=randn(T,1);                       % GPi Amplitude(?)
y=randn(T,1);                       % STN Amplitude(?)
z=angle(exp(1i*2*pi*rand(T,1)));    % Phase-difference(?)

% Truncate data range to exclude sparse regions - amplitude only
stdrange = 2;
if (~isempty(stdrange))
    % Mean subtract
    x = x-mean(x);
    y = y-mean(x);
    % Isolate deviant samples and reject
    nokeep = (abs(x)>stdrange*std(x)) | (abs(y)>stdrange*std(y));
    x = x(~nokeep);
    y = y(~nokeep);
    z = z(~nokeep);
    T = length(x);
end;

% Stratify --- uniform across range, so transform first if doing so (i.e. log-transform, etc.)
levels = [ 10 10 20 ];      % <-- a smaller stratification range can help sparse sampling
% Histogram without bias on x (amplitude)
x0 = x; x = ones(size(x0));     % use `ones` so that (x0>xq(k)) does not have to be >= which double counts
xq = quantile(x0,(0:1/levels(1):1.0));
for k = (1:levels(1))
    x((x0>xq(k)) & (x0<=xq(k+1))) = k;
end;
% Histogram without bias on y (amplitude)
y0 = y; y = ones(size(y0));
yq = quantile(y0,(0:1/levels(1):1.0));
for k = (1:levels(1))
    y((y0>yq(k)) & (y0<=yq(k+1))) = k;
end;
% Histogram without bias on z (using phase wrapping)
z=round((z-min(z))/range(z)*levels(3)+1); z(z==(levels(3)+1)) = 1;
phases = (-pi:2*pi/levels(3):(pi-2*pi/levels(3)))+pi/levels(3);

% Construct joint distribution (per sample)
pxyz=zeros(levels);
for t=(1:T)
    pxyz(x(t),y(t),z(t)) = pxyz(x(t),y(t),z(t)) + 1;
end;
% 3D convolution -- this is a way of blurring the sparse sampling issue
%  warning: introduces edge effects in its current implementation
if ( false )
    pxyz = convn( pxyz, ones(5), 'same' );
end;
% Histogram -> prob dist
pxyz = pxyz/sum(pxyz(:));

% MI per phase-difference bin
I = zeros(levels(3),1);
for n = (1:levels(3))
    I(n) = mi( pxyz(:,:,n) );
end;

% TE per phase-difference bin
TExy = zeros(levels(3),1);
TEyx = TExy;
for n = (1:levels(3))
    TExy(n) = te( y(z==n), x(z==n) );
    TEyx(n) = te( x(z==n), y(z==n) );
end;

% MI conditioned on phase-difference (z)
Icond = mi_cond( pxyz );

% Plot
figure;
precision = 1e6;
subplot(221);
    cla; plot( phases, I ); hold on;
    ylabel('MI(x;y)'); xlabel('\Delta\Phi');
    h=plot(mean(xlim),mean(ylim),mean(xlim),mean(ylim),mean(xlim),mean(ylim),'linestyle','none');
    title(['Conditional mutual information MI(x;y|\Delta\Phi): ' num2str(round(precision*Icond)/precision)]);
subplot(222);
    cla; h = plot( phases, TExy, phases, TEyx ); hold on;
    ylabel('TE_{a\rightarrowb}'); xlabel('\Delta\Phi');
    legend(['TE_{x\rightarrowy} (' num2str(round(precision*te(y,x))/precision) ')'],['TE_{x\rightarrowy} (' num2str(round(precision*te(x,y))/precision) ')']);
    title(['Transfer entropies']);
subplot(224);
    cla; plot( phases, TExy-TEyx ); hold on; plot( xlim, [0 0], 'k' );
    ylim(max(abs(ylim))*[-1 1]); xlabel('\Delta\Phi');
    ylabel('TE_{x\rightarrowy} - TE_{y\rightarrowx}');
    legend(['x\rightarrowy > y\rightarrowx (' num2str(round(precision*(te(y,x)-te(x,y)))/precision) ')']);
    title('Transfer entropy directionality');
