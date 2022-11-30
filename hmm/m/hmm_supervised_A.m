function A = hmm_supervised_A(supervised_states,statecount)
%A = hmm_supervised_A(supervised_states,statecount)

% Index from 1 not 0
if (min(supervised_states)==0)
    supervised_states = supervised_states + 1;
end;
if (~exist('statecount'))
    statecount = max(supervised_states);
end;

A=zeros(statecount);
for n1=(1:statecount)
    for n2=(1:statecount)
        A(n1,n2) = length(find( (supervised_states(1:(end-1))==n1) & (supervised_states(2:end)==n2) ));
    end;
end;
A=hmm_rownorm(A);
