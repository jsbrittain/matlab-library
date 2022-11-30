function [z,x,x1,x2] = modelS2fit( attained, target, clamptrials, permutations )

% Optimisation options
options = optimset(@fmincon);
options = optimset(options,'TolCon',1e-6);
options = optimset(options,'Algorithm','interior-point');

% Model constraints
x0 = [0.6 0.2 0.9924 0.02];
A = [ 1  0 -1 0 ;       % Af < As
      0 -1  0 1 ];      % Bf > Bs
b = [ 0; 0 ];
Aeq = []; beq = []; nlcon = [];
LB = zeros(size(x0)); UB = ones(size(x0));

% Precompute mean and stddev
if (size(attained,1)>1)
    switch ( 2 )
        case 1,     % Sum-of-squares
            attainedvar =  ones(1,size(attained,2));
        case 2,     % Chi-square goodness-of-fit
            attainedvar =  nanstd( attained, [], 1 ).^2;
    end;
else
    attainedvar = ones(size(attained));
end;
attained   = nanmean( attained, 1 );
target     = nanmean( target,   1 );

% Iterate permutations
for n = (1:permutations)
    x0i = rand(size(x0));
    try
        [z(:,n),fval,exitflag,output] = fmincon( @(x)modelS2err(x,attained,attainedvar,target,clamptrials), x0i, A, b, Aeq, beq, LB, UB, nlcon, options );  % Constrained
    catch
        z(:,n) = NaN;
    end;
end;

% Freerun model
[x,x1,x2] = modelS2( nanmedian(z,2), target, clamptrials );
