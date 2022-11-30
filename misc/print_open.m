function print_open
% function print_open
%
% Prints all open figures to the default printer
%

fh=findobj('type', 'figure');
for ind=1:length(fh)
    print(fh(ind));
end;
