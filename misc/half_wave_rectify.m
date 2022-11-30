function x = half_wave_rectify(x);

ind=find(x<0);
x(ind)=0;
