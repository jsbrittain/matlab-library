function triallist = epoch_trim_dual( erptime, epochs1, epochs2, strTitle )

meandisplay = true;

% Launch and shift focus to 
fh = figure;
set(gcf,'windowstyle','docked');
figure(fh);
plotedit('on');

% Axes
ah(1) = subplot(1,2,1); hold('on'); axis(ah(1),'off'); title(strTitle);
ah(2) = subplot(1,2,2); hold('on'); axis(ah(2),'off');
set(fh,'color',[1 1 1]);

% Plot epochs individually
axes(ah(1));
for n = (1:size(epochs1,3))
    plot(erptime,epochs1(:,1,n),'color',[0 0 1],'linewidth', 1);
end;
axes(ah(2));
for n = (1:size(epochs1,3))
    plot(erptime,epochs2(:,1,n),'color',[0 0 1],'linewidth', 1);
end;

% Record all trace handles
h1all = get(ah(1),'children');
h2all = get(ah(2),'children');

h1mean = zeros(0,3);
h2mean = zeros(0,3);

% Wait for signal to end
str = '';
while ( ~strcmpi(str,'x') )
    % Update mean
    axes(ah(1)); delete(h1mean);
    h1mean(1) = plot( erptime, mean(cell2mat(get(get(ah(1),'children'),'ydata')),1), 'r-' );
    h1mean(2) = plot( erptime, mean(cell2mat(get(get(ah(1),'children'),'ydata')),1)-2*std(cell2mat(get(get(ah(1),'children'),'ydata')),[],1), 'r-' );
    h1mean(3) = plot( erptime, mean(cell2mat(get(get(ah(1),'children'),'ydata')),1)+2*std(cell2mat(get(get(ah(1),'children'),'ydata')),[],1), 'r-' );
    axes(ah(2)); delete(h2mean);
    h2mean(1) = plot( erptime, mean(cell2mat(get(get(ah(2),'children'),'ydata')),1), 'r-' );
    h2mean(2) = plot( erptime, mean(cell2mat(get(get(ah(2),'children'),'ydata')),1)-2*std(cell2mat(get(get(ah(2),'children'),'ydata')),[],1), 'r-' );
    h2mean(3) = plot( erptime, mean(cell2mat(get(get(ah(2),'children'),'ydata')),1)+2*std(cell2mat(get(get(ah(2),'children'),'ydata')),[],1), 'r-' );
    % Query user for input
    str = input(' [R]emove trace, [H] to toggle highlight, e[X]it ... ', 's');
    % Check if input makes sense
    if (strcmpi(str,'h') || strcmpi(str,'r'))
        % Find object on appropriate axis
        h_index = find( h1all == gco, 1, 'first' );
        h_index = [ h_index find( h2all == gco, 1, 'first' ) ];
        % Highlight selection on both axes
        if (strcmpi(str,'h'))
            if (isequal(get(gco,'color'),[1 0 0]))
                set( h1all(h_index), 'color', [0 0 1], 'linewidth', 1 );
                set( h2all(h_index), 'color', [0 0 1], 'linewidth', 1 );
            else
                set( h1all(h_index), 'color', [1 0 0], 'linewidth', 2 );
                set( h2all(h_index), 'color', [1 0 0], 'linewidth', 2 );
            end;
        end;
        % Delete trials
        if (strcmpi(str,'r'))
            delete( h1all(h_index) );
            delete( h2all(h_index) );
        end;
    end;
end;
delete(h1mean);
delete(h2mean);

% Compare remaining trace handles to initial set
h1remaining = get(ah(1),'children');
triallist = flipud( ismember( h1all, h1remaining ) );

% Close figure
close( fh );
