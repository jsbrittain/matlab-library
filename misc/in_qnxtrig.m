function [trig] = in_qnxtrig(filestem)
%  function [trig] = in_qnxtrig(filestem)
%
%  Reads in trigger data from QNX generated trigger file
%  Arguments - filestem  : filename (stem)
%  Returns   - trig      : Trigger data in columns
%              wf_fil    : W filestem
%              wf_no     : W file number 

% This version 1.00, 25/11/04, DMH.

if (nargin<1)
  error(' Not enough input arguments');
end  

if (nargout<1)
  error(' Not enough output arguments');
end  

% Open QNX trigger file
filname=sprintf('%s.txt',filestem);
fid=fopen(filname,'r');
if (fid==-1)
   error(['Cannot open file: ',filname]);
end
disp(['Reading: ',filname]);

% Read and disply header line, ATF version no.
run_str=fscanf(fid,'%s',1);
wf_fil=fscanf(fid,'%s',1);
wf_str=fscanf(fid,'%s',1);
wf_no=fscanf(fid,'%d',1);
line=fgetl(fid);                    % Read in to end of line
line1=fgetl(fid);                   % Read in header line 1
line2=fgetl(fid);                   % Read in header line 2

% Read in trigger and conver to column format
trig=fscanf(fid,'%f,%f',[2,inf])'; 

fclose(fid);

disp(['Data File: ',wf_fil,'  W file no: ',num2str(wf_no),',  Triggers: ',num2str(length(trig))]);
