function [X,Y,Z] = rotMatricesDeg( rotx, roty, rotz )

rotx = rotx * 2 * pi / 360;
roty = roty * 2 * pi / 360;
rotz = rotz * 2 * pi / 360;

[X,Y,Z] = rotMatrices( rotx, roty, rotz );
