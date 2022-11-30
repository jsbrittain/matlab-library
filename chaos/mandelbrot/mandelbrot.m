%% Open multicore processing pool

matlabpool open

%% Compute Mandelbrot set
%
% z -> z^2 + c
%

tic

N = 1000; D = 4;

xs = (-2:0.01:0.4);
ys = (-1:0.01:1);

M = zeros(length(xs),length(ys));

parfor xn = (1:length(xs))
    
    Mm = M(xn,:);
    
    z = zeros(1,N);
    for yn = (1:length(ys))
        
        c = xs(xn) + 1i*ys(yn);
        
        n = 2;
        while ((n<=N) && (abs(z(n-1))<D))
            z(n) = z(n-1)^2 + c;
            n = n + 1;
        end;
        
        Mm(yn) = (n - 1);
    end;
    
    M(xn,:) = Mm;
end;

imagesc(xs,ys,log(M'));

toc

matlabpool close