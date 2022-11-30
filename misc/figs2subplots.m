function newfig = figs2subplots(handles,arr);
% FIGS2SUBLPLOTS Combine axes in many figures into subplots in one figure
%
%   The syntax:
%
%       >> newfig = figs2subplots(handles,arr);
%   
%   creates a new figure with handle "newfig", in which the axes specified
%   in vector "handles" are reproduced and aggregated as subplots. 
%
%   Vector "handles" is a vector of figure and/or axes handles. If an axes
%   handle is encountered, the corresponding axes is simply reproduced as
%   a subplot in the new figure; if a figure handle is encountered, all its
%   children axes are reproduced as subplots in the figure.
%
%   Vector "arr" is an optional subplot arrangement vector of the form 
%   [M N], where M and N specify the number of rows and columns of
%   subplots to use. When the k-th handle is processed, the command
%   "subplot(M,N,k)" is issued, unless k > MN (in which case the current
%   handle and all those following it are ignored). Default for "arr" is
%   [Na 1], where Na is the total number of axes found while parsing vector
%   "handles".

% Parsing handles vector
av = [];
for k = 1:length(handles)
    if strcmp(get(handles(k),'Type'),'axes')
        av = [av handles(k)];
    elseif strcmp(get(handles(k),'Type'),'figure');
        fc = get(handles(k),'Children');
        for j = length(fc):-1:1
            if strcmp(get(fc(j),'Type'),'axes')
                av = [av fc(j)];
            end;
        end;
    end;
end;

% Setting the subplots arrangement
Na = length(av);
if nargin < 2
    arr = [Na 1];
    Ns = Na;
else
    Ns = prod(arr);
end;

% Creating new figure
da = zeros(1,Ns);
newfig = figure;
for k = 1:min(Ns,Na)
    da(k) = subplot(arr(1),arr(2),k);
    na = copyobj(av(k),newfig);
    set(na,'Position',get(da(k),'Position'));
    delete(da(k));
end;