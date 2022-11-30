function epochs=pnormepochs(epochs,chlist)

for ind=(1:size(epochs,3))
    for ch=chlist
        epochs(:,ch,ind)=epochs(:,ch,ind)/sqrt(mean(abs(epochs(:,ch,ind)).^2));
    end;
end;
