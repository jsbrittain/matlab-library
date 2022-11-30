function out=downsample_jsb(dat,oldrate,newrate)
%function out=downsample_jsb(dat,oldrate,newrate)
%
% Downsample
%
%function out=downsample_jsb(dat,oldrate,newrate)

if (oldrate==newrate)
    out=dat;
    return;
end;

% Down-sample signals
windowSize=oldrate/newrate;
[b,a]=butter(4,(newrate/2.5)/(oldrate/2),'low');
dat=filtfilt(b,a,dat);
samples=size(dat,1);
k=round(windowSize:windowSize:samples);
out=dat(k,:);
