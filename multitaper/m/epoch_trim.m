function triallist = epoch_trim( erptime, epochs )
%function triallist = epoch_trim( erptime, epochs )
%
%function triallist = epoch_trim( erptime, epochs )

% Launch and shift focus to 
fh = figure;
set(gcf,'windowstyle','docked');
figure(fh);

% Axis
h = subplot(1,1,1);
hold('on');

% Plot epochs individually
for n = (1:size(epochs,3))
    plot(erptime,epochs(:,1,n),'color',[0 0 1],'linewidth', 1);
end;

% Record all trace handles
h0 = get(h,'children');
axis(h,'off');
set(fh,'color',[1 1 1]);

% Wait for signal to end
str = ' ';
while ( ~isempty(str) )
    str = input(' Remove unwanted traces ([H] to toggle highlight), then hit [ENTER]...', 's');
    if (strcmpi(str,'h'))
        if (isequal(get(gco,'color'),[1 0 0]))
            set( gco, 'color', [0 0 1], 'linewidth', 1 );
        else
            set( gco, 'color', [1 0 0], 'linewidth', 2 );
        end;
    end;
end;

% Compare remaining trace handles to initial set
h1 = get(h,'children');
triallist = flipud( ismember( h0, h1 ) );

% Close figure
close( fh );
