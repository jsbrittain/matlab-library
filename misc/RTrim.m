function out = RTrim(in_str)
% function out = RTrim(in_str)
%
% Trim trailing spaces from a string
%

last_point=0;
for ind=length(in_str):-1:1
    if in_str(ind)~=' '
        last_point = ind;
        break;
    end;
end;

if (last_point == 0)
    out = '';
else
    out(1:last_point) = in_str(1:last_point);
end;
