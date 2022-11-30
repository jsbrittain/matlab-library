function sp = sp_baseline( sp, spect )

if (size(spect,2)==1)
    spect = repmat( spect, 1, size(sp,3) );
end;

for n = (1:size(sp,3))
    sp(:,:,n) = sp(:,:,n) ./ repmat(spect(:,n),1,size(sp,2));
end;
