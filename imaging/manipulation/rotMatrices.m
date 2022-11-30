function [X,Y,Z] = rotMatrices( rotx, roty, rotz )

X = [ 1 0 0 0;
      0 cos(rotx) sin(rotx) 0
      0 -sin(rotx) cos(rotx) 0
      0 0 0 1 ];
  
Y = [ cos(roty) 0 -sin(roty) 0;
      0 1 0 0;
      sin(roty) 0 cos(roty) 0
      0 0 0 1 ];
  
Z = [ cos(rotz) sin(rotz) 0 0;
      sin(rotz) cos(rotz) 0 0;
      0 0 1 0;
      0 0 0 1 ];
