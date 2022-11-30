function figsave_open(filstem)
% Save the group of open figures to filstem<no.>

% Get figure handles
figs=sort(findobj(0,'type','figure'))';

% Recurse figures and save
for ind=1:length(figs)
    saveas(figs(ind),[filstem int2str(ind)]);
end;
