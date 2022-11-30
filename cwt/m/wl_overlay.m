function wl_overlay( freqs, time, p, alpha, flip, plotmethod, border )
%function wl_overlay( time, freqs, p, alpha, flip, plotmethod, border )
%
% for use with wl_ttest, wl_ttest2
%
%       flip                should be 'true' for use with wlpsp routines
%       border              add border (p=1) to image, required for correct
%                            derivation of contours in some cases [true]
%       plotmethod          1 - outline, 2 - opaque, 3 - translucent
%
%function wl_overlay( time, freqs, p, alpha, flip, plotmethod, border )

% Check input arguments
if (~exist('alpha','var'))
    alpha = [];
end;
if (~exist('flip','var'))
    flip = [];
end;
if (~exist('plotmethod','var'))
    plotmethod = [];
end;
if (~exist('border','var'))
    border = [];
end;

% Default parameters
if (isempty(alpha))
    alpha = 0.05;
end;
if (isempty(flip))
    flip = false;
end;
if (isempty(plotmethod))
    plotmethod = 0;
end;
if (isempty(border))
    border = true;
end;
if (isempty(time))      % Needs to be after `flip'
    warning(' Using default time axis!');
    time = ylim;
    time = time(1):(time(end)-time(1))/(size(p,1)-1):time(end);
end;
if (isempty(freqs))      % Needs to be after `flip'
    warning(' Using default frequency axis!');
    freqs = [0 1];%xlim;
    freqs = freqs(1):(freqs(end)-freqs(1))/(size(p,2)-1):freqs(end);
end;

% Flip axes (provides compatibility with wlpsp routines)
if ( flip )
    p = p';
end;

% Add border of insignificance around p-matrix
if ( border )
    p = [ ones(1,size(p,2)+2); [ ones(size(p,1),1) p ones(size(p,1),1) ]; ones(1,size(p,2)+2)  ];
    dt=time(2)-time(1); df=freqs(2)-freqs(1);
    freqs = [freqs(1)-df freqs freqs(end)+df];
    time = [time(1)-dt time time(end)+dt];
end;

% Overlay
hold('on');
switch (plotmethod)
    
    case 0,     % Outline
        contour(freqs,time,(p<alpha),1,'k');
        
    case 1,     % Opaque
        % Get internal contours
        c=contourc(freqs,time,double(p<alpha),1);
        k = 1;
        while (k<size(c,2))
            n = c(2,k);
            patch(c(1,k+(1:n)),c(2,k+(1:n)),1,'facecolor',[1 1 1]);
            k = k + n + 1;
        end;
        
    case {2,3,4},     % Translucent
        % Get internal contours
        c=contourc(freqs,time,double(p<alpha),1);
        cx = {}; cy = {};
        m = 1; k = 1;
        while (k<size(c,2))
            n = c(2,k);
            cx{m} = c(1,k+(1:n));
            cy{m} = c(2,k+(1:n));
            k = k + n + 1;
            m = m + 1;
        end;
        % Determine patch objects
        if ( false )
            % Compute contour patches
            [f, v] = poly2fv(cx,cy);
        else
            % Outer boundary (this effectively inverts the filled region)
            xlims = xlim; ylims = ylim;
            x0 = xlims([0 1 1 0]+1); y0 = ylims([0 0 1 1]+1);
            % Compute inverted contour patches
            [f, v] = poly2fv({x0,cx{:}},{y0,cy{:}});
        end;
        % Display patches
    	switch (plotmethod)
            case 2, patch( 'Faces', f, 'Vertices', v, 'FaceColor', [1 1 1], 'EdgeColor', 'none' );
            case 3, patch( 'Faces', f, 'Vertices', v, 'FaceColor', [1 1 1], 'EdgeColor', 'none', 'facealpha', 0.5 );
            case 4, patch( 'Faces', f, 'Vertices', v, 'FaceColor', [1 1 1], 'EdgeColor', 'none', 'facealpha', 0.5 );
                    contour( freqs,time,double(p<alpha),1,'linecolor',[0 0 0],'linewidth',2 );
                    
        end;
        
    otherwise
        error(' Unknown plotmethod');
        
end;
