function [states,reorder] = hmm_qfit(oldstates,q)

Q = max([oldstates q]);
states = oldstates;

% Determine best fit vectors
len = zeros(Q);
for n1=(1:Q)
    for n2=(1:Q)
        len(n1,n2) = length(find(q(oldstates==n1)==n2));        
    end;
    [dummy,reorder]=max(len,[],2);
end;

% Re-order state vector
for n1=(1:Q)
    states(oldstates==n1)=reorder(n1);
end;
