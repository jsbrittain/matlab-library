function p=pdist_pref(epsilon,X,Xi,jj,w1,w2,k,Pref)
%function p=pdist_pref(epsilon,X,Xi,jj,w1,w2,k,Pref)
%
% Error function for epsilon search
% For use in pdist_* functions
%
% Minimises
%   | (Stam 2005, Eq.2) - Pref |
%
%function p=pdist_pref(epsilon,X,Xi,jj,w1,w2,k,Pref)

p=abs(sum(((epsilon-sqrt(sum((Xi-X(jj,:,k)).^2,2)))>0),1)/(2*(w2-w1))-Pref);
