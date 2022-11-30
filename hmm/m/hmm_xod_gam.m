function [Dpdf,Dparams,Dmax,Drnd] = hmm_xod_gam(state,mindist)
%function [Dpdf,Dparams,Dmax,Drnd] = hmm_xod_gam(state)

% Determine input arguments
states = max(state);
if (~exist('mindist'))
    mindist=[];
end;
if (isempty(mindist))
    mindist=Inf;
end;

% Explicit state duration statement
for k=(1:states)
    state_start = find(diff(state==k)>0);
    state_end = find(diff(state==k)<0);
    state_start = state_start( state_start < state_end(end) );
    state_end = state_end( state_end > state_start(1) );
    state_dur = (state_end - state_start);
    state_dur(state_dur<mindist) = [];
    
    Dpdf{k} = @gampdf;
    Drnd{k} = @gamrnd;
    Dparams{k} = num2cell(gamfit( state_dur ));
    Dmax(k) = find((cumsum(gampdf([1:1:10000],Dparams{k}{:})) + gampdf(1,Dparams{k}{:}) ) > 0.99,1,'first');
end;
