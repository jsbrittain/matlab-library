function c = ann_cost_vec( xInput, wvec, bvec, g, y, nNeurons )

w = cell(1,length(nNeurons)-1);
b = cell(1,length(nNeurons)-1);

kw = 1;
kb = 1;

for layerNo = (1:(length(nNeurons)-1))
    wn = [ nNeurons(layerNo) nNeurons(layerNo+1) ];
    w{layerNo} = reshape( wvec(kw:(kw+wn(2)*wn(1)-1)), wn(2), wn(1) );
    kw=kw+wn(2)*wn(1);
    b{layerNo} = reshape( bvec(kb:(kb+wn(2)-1)), wn(2), 1 );
    kb=kb+wn(2);
end

c = ann_cost( xInput, w, b, g, y );
