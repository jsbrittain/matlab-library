function [events,triggers]=get_triggers(trig,rate)
%function [events,triggers]=get_triggers(trig,rate)
%
% Input parameters
%	trig	Transition lines (low-to-high bits)
%	rate	Sampling rate
%
%function [events,triggers]=get_triggers(trig,rate)

% Preprocess (nb: some channels may have no events)
M=size(trig,2);                 % Trigger channels

trig=trig-min(trig(:));
trig=trig/max(trig(:));
trig1=zeros(size(trig));
trig1(trig>0.5)=1;
trig=trig1;

trig=diff(trig,[],1);			% Differentiate trig channels
trig=trig/max(trig(:));			% Normalise trigger diffs

% Isolate low-high transitions
for ind=(1:M)
    trans{ind}=find(trig(:,ind)>0.5);
end;

% Construct trigger channel (permitting 5msec error between channels)
triggers=zeros(size(trig,1),1);
triggers(trans{1})=1;
% Add trigger lines (summate if pulse exists within 5ms)
for ch=(2:M)
    bit=2^(ch-1);
    triglist=find(triggers);
    for ind=(1:length(trans{ch}))
	index=find(abs(triglist-trans{ch}(ind))<=(0.002*rate));
	if (~isempty(index))
	    triggers(triglist(index))=triggers(triglist(index))+bit;
	else
	    triggers(trans{ch}(ind))=bit;
	end;
    end;
end;

% Generate events structure
for ind=(1:((2^M)-1))
    events{ind}=find(triggers==ind);
end;
