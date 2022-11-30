function eventtrigs = epoch_peruse(dat,eventtrigs,rate,offset,duration, baseline, markall)
%
%
% Options as epochdata
%   NB: Takes a SINGLE data channel only
%
% Last option: markall, mark all trials, disables checking
%

if (~exist('markall','var'))
    markall = [];
end;
if (isempty(markall))
    markall = false;
end;

[epochs,erptime] = epochdata(dat,eventtrigs,rate,offset,duration, baseline);

zeropt = dsearchn(erptime',0);
rejecttrials = [];
channel = 1;

fh = figure;
k = 1;
while (k <= length(eventtrigs))
    % Plot trial(s)
    cla; hold('on');
    h = plot(squeeze(epochs(:,channel,:))); set(h,'color',[1 1 1]*0.8);
    for n = (1:length(rejecttrials))
        h = plot(epochs(:,channel,rejecttrials(n))); set(h,'color',[1 0.4 0.4]);
    end;
    if (ismember(k,rejecttrials))
        plotcol = 'r';
    else
        plotcol = 'b';
    end;
    plot( epochs(:,channel,k), plotcol );
    plot(xlim,[0 0],'k'); plot(zeropt*[1 1],ylim,'k');
    title(['TRIAL ' num2str(k) ' of ' num2str(length(eventtrigs))]);
    % Query next action
    if (markall)
        reply = 'C';
    else
        reply = input('(C)orrect, (B)ack, (D)elete/undelete, e(X)it, or accept [ENTER]: ', 's');
    end;
    if (~isempty(reply))
        if (strcmpi(reply,'C'))
            % Get new trigger and replot
            [x,y] = ginput(1);
            if (isempty(y))
                continue;
            end;
            eventtrigs(k) = eventtrigs(k) + round(x) - zeropt;
            [epochs,erptime] = epochdata(dat,eventtrigs,rate,offset,duration, baseline);
            if (markall)
                k = k + 1;
            end;
        elseif (strcmpi(reply,'B'))
            % Go back one
            if (k>1)
                k = k - 1;
            end;
        elseif (strcmpi(reply,'D'))
            % Delete / undelete trial
            if (ismember(k,rejecttrials))
                rejecttrials = setdiff(rejecttrials,k);
            else
                rejecttrials = unique([rejecttrials k]);
            end;
        elseif (strcmpi(reply,'X'))
            % Abort
            break;
        else
            disp(' Unknown option!');
        end;
        % Decrement trial count (keeps the same)
        k = k - 1;
    else
        if (k==length(eventtrigs))
            msgbox('Finished, closing window.','Finished');
        end;
    end;
    % Next trial
    k = k + 1;
end;
close(fh);

eventtrigs = eventtrigs( setdiff(1:length(eventtrigs),rejecttrials) );
