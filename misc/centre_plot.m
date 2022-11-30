function centre_plot(ah);
%function centre_plot(ah);
%
% Centralise axis (with colorbar if one exists)
%
% Input parameters
%       ah      Axis handle (default: gca)
%
%function centre_plot(ah);

% Check input parameters
if (nargin<1)
    ah=gca;
end;
if (~ishandle(ah))
    error(' Handle provided does not point to a valid graphic object.');
end;

% Determine if a colorbar is associated with the axis
ch=[];      % Colorbar handle
cbhandles=get(findall(gcf,'type','image','tag','TMW_COLORBAR'),{'parent'});
for ind=1:length(cbhandles)
    ud=get(cbhandles{ind},'userdata');
    if (isfield(ud,'PlotHandle'))
        if (ud.PlotHandle==ah)
            ch=cbhandles{ind};
            break;
        end;
    end;
end;

% Get axis position
ahpos=get(ah,'position');
if (isempty(ch))
    chpos=[0 0 0 0];
else
    chpos=get(ch,'position');
end;

% Reposition axis
width=ahpos(3)+max(0,chpos(1)+chpos(3)-ahpos(1)-ahpos(3));
shift=0.5-width/2-ahpos(1);
set(ah,'position',[ahpos(1)+shift ahpos(2:4)]);

% Reposition colorbar if it exists
if (~isempty(ch))
    set(ch,'position',[chpos(1)+shift chpos(2:4)]);
end;
