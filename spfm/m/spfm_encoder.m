function v=spfm_encoder(v,x,tau,dt);
%function v=spfm_encoder(v,x,tau,dt);
%
% Sigma Pulse Frequency Modulation (SPFM) Encoder
%
% Discrete implementation of the encoder using a backward Euler integration
% scheme (as Ref).
%
% Input parameters
%       v       Output of encoder at time index n
%       x       Input to encoder at time index n
%       tau     Time constant
%       dt      Sampling interval
%
% Output parameter
%       v       Output for time index (n+1)
%
% Ref: D.M. Halliday (1998) "Generation and characterization of correlated
%       spike trains". Computers in Biology and Medicine (28) pp.143-152
%
%function v=spfm_encoder(v,x,tau,dt);

v=(v+(dt/tau)*x)/(1+dt/tau);
