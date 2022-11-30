function dpsssave(N,W)
% A simple utility function to calculate dpss and save them to disk
% with filename and format as used by the various mtm routines.
% 
% Alternatively, a script something like the following could be
% used to totally automate the process - run this one overnight
% if you are going to compute many large dpss!
%
% cd e:\matlab\mtm    % Substitute your desired directory here.
% [E,V]=dpsscalc(256,2);   save s256_2.mat
% [E,V]=dpsscalc(256,5/2); save s256_25.mat
% [E,V]=dpsscalc(256,3);   save s256_3.mat
% [E,V]=dpsscalc(256,7/2); save s256_35.mat
% [E,V]=dpsscalc(512,2);   save s512_2.mat
% [E,V]=dpsscalc(512,5/2); save s512_25.mat
% [E,V]=dpsscalc(512,3);   save s512_3.mat
% [E,V]=dpsscalc(512,7/2); save s512_35.mat
%
% etc., etc.
%
% Written by Eric Breitenberger, version date 10/1/95.
% Please send comments and suggestions to eric@gi.alaska.edu
%


dir=input('Directory where DPSS are to be stored: ', 's');
eval(['cd ' dir]);

if rem(W,1)==0
  filename=['s' int2str(N) '_' int2str(W) '.mat']; 
else
  filename=['s' int2str(N) '_' int2str(10*W) '.mat']; 
end

disp('Calculating the DPSS...')

eval(['save ' filename])


