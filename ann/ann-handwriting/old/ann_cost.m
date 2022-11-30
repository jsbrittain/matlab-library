function c = ann_cost( x, w, b, g, y )

for layerNo = (1:(length(w)))
    x{layerNo+1} = g( w{layerNo}*x{layerNo} + b{layerNo} );
end
c = mean( mean((y - x{end}).^2) );
