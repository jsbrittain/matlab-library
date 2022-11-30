function [t,pval1,pval2] = analysisDriftPermPlot( metric, nulldistr, channames )

chcount = size(nulldistr,1);
permutations = size(nulldistr,2);

% p=value
for ch = (1:chcount)
    p1 = 1-find((metric(ch)<=sort(nulldistr(ch,:))),1,'first')/permutations;        % One-sided (assumes metric +ve)
    if (isempty(p1)), p1 = 1/permutations; end; if (p1>0.5), p1=1-p1; end;
    pval1(ch) = p1;
    p2 = [ 1-find((abs(metric(ch))<=sort(nulldistr(ch,:))),1,'first')/permutations + find((-abs(metric(ch))<=sort(nulldistr(ch,:))),1,'first')/permutations ];        % Two-sided
    if (isempty(p2)), p2 = 1/permutations; end;
    pval2(ch) = p2;
end;

% Plot permutation distributions
figure;
rows = ceil(sqrt(chcount));
cols = ceil(chcount/rows);
t = nan(chcount,1);
for ch = (1:chcount)
    subplot(rows,cols,ch);
        [count,bins] = hist( nulldistr(ch,:) );
        count = count/sum(count)/mean(diff(bins));
        bar( bins, count ); hold('on');
        mu = mean(nulldistr(ch,:));
        sd = std(nulldistr(ch,:));
        t(ch) = (metric(ch)-mu)/sd;
        plot( bins, normpdf(bins,mu,sd), 'c--' );
        plot( metric(ch)*[1 1], ylim, 'r', 'linewidth', 2 );
        xlabel(sprintf('t = %4.3f, p = %4.3f (p = %4.3f one-sided)',t(ch),pval2(ch),pval1(ch)));
        title(channames{ch});
end;
