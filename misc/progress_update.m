function progress = progress_update( progress, k, N )

newprogress = round(100*k/N);
if ( newprogress > progress )
    progress = newprogress;
    if ((mod(progress,10)==0) && (progress~=100))
        fprintf('%g',progress);
    else
        fprintf('.');
    end;
end;
