function xlimits(arg1,arg2,arg3);
%function xlimits([fh],start,stop);
%
% Recurse subplots of current figure and change x-range
%
% Input parameters
%       fh          (opt) Figure handle
%       start       Minimum x-limit
%       stop        Maximum x-limit
%
%function xlimits([fh],start,stop);

% Determine input parameters
switch (nargin)
    case 2, fh=gcf;  start=arg1; stop=arg2;
    case 3, fh=arg1; start=arg2; stop=arg2;
    otherwise
        error('Incorrect input parameters.');
end;

% Determine subplot vector (exclude colorbars)
ah=get(fh,'children');
ch=get(findall(fh,'type','image','tag','TMW_COLORBAR'),{'parent'});
ah=setdiff(ah,[ch{:}]);

% Recurse subplots and change x-limits
for ind=1:length(ah)
    subplot(ah(ind));
    xlim([start stop]);
end;
