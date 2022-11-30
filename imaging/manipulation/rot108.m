% Translate standard space MNI152 , rotate a set amount and recentre

T = transMatrix( -90, -126, -72 );
[X,Y,Z] = rotMatricesDeg( 19, 0, 0 );

% Transpose ALL matrices into FSL convention
T=T'; X=X'; Y=Y'; Z=Z';

% Form composite matrix
M = inv(T)*X*Y*Z*T;

% Inverse transform
%M = inv(M);

% Write to text file
fid=fopen('rotmat.mat','w');
fprintf(fid, '%g %g %g %g\n%g %g %g %g\n%g %g %g %g\n%g %g %g %g\n',M');
fclose(fid);
