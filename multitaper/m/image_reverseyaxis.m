function image_reverseyaxis(ah)
%function image_reverseyaxis([ah])
%
% Reverse the tick labels on the y axis for plots created
% with the 'image' function
%
% Input parameter
%   ah      (opt) Axis handle
%

if (nargin==0)
    ah=gca;
end;

oldytick=get(ah,'ytick');
oldyticklabel=get(ah,'yticklabel');
ylims=ylim;
newytick=fliplr((ylims(2)-oldytick+ylims(1)));
set(ah,'ytick',newytick);
set(ah,'yticklabel',flipud(oldyticklabel));
