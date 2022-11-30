function stats = diffOfCorrelationTest( ra, rb, na, nb )
% Note: This test is applied to the sample correlation coefficient "r",
% NOT the coefficient of determination R^2!

za = atanh( ra );
zb = atanh( rb );

vara = na/(na-3);
varb = nb/(nb-3);

stats.z = ( za - zb )/sqrt( vara/na + varb/nb );
stats.pRight = normcdf( stats.z );
stats.pLeft = 1-normcdf( stats.z );
stats.pTwoTailed = 2*(1-normcdf( abs(stats.z) ));
