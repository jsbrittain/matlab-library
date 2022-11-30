function c = ann_cost_vec_optim( vec, xInput, g, y, nNeurons )

bcount = sum(nNeurons(2:end));
wvec = vec(1:(length(vec)-bcount+1));
bvec = vec((length(vec)-bcount+1):end);

c = ann_cost_vec( xInput, wvec, bvec, g, y, nNeurons );
