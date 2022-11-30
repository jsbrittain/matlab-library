%function whatcha = vars2struct

vnames = who;
whatcha = struct;

for n=(1:length(vnames))
    whatcha.(vnames{n}) = eval(vnames{n});
end;
