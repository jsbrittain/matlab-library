function smr_stackplot(dat,channels)



figure;
stackplot(time,single([dat{channels}]));
