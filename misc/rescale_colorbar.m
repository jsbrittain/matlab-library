function ch=rescale_colorbar(ah);
%function rescale_colorbar(ah);
%
% Reduce influence of colorbar on axis.
% Creates a new colobar if one is not currently associated with the axis.
%
% Input parameters
%       ah      master axis handle (default: gca)
%
%function rescale_colorbar(ah,ch);

% Check input parameters
if (nargin<1)
    ah=gca;
end;
if (~ishandle(ah))
    error(' Handle provided does not point to a valid graphic object.');
end;

% Determine if colorbar exists, otherwise create one (preserves existing formating)
ch=[];
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
if (isempty(ch))
    ch=colorbar('peer',ah);
end;

% Rescale axes to reduce colorbar influence
ahpos=get(ah,'position');
chpos=get(ch,'position');
width=1/4;                      % Rescale to a portion of original colorbar width
ahpos(3)=ahpos(3)+(1-width)*chpos(3)+(chpos(1)-ahpos(1)-ahpos(3))/2;
chpos(1)=chpos(1)+(1-width)*chpos(3);
chpos(3)=width*chpos(3);
set(ah,'position',ahpos);
set(ch,'position',chpos);
