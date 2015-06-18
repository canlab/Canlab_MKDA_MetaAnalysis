function [u,xpdf,ypdf,xcdf,ycdf,h] = d_pdf(vec,crit_p,varargin)
% [u,xpdf,ypdf,xcdf,ycdf,h] = d_pdf(vec,crit_p,[opt] plot, [opt] figure handles)
%
% vec       vector of values with unknown distribution
% crit_p    p value at which to make height threshold (e.g., .05)
% plot      Optional.  0 = none, 1 = pdf, 2 = pdf + cdf
%           followed by plot handles
% h         handles to line plots on each figure
%
% Tor Wager, 2/6/02

doplot = 0;
if length(varargin) > 0
    doplot = varargin{1};
    for i = 2:length(varargin)
        fh(i-1) = varargin{i};
    end
end

% -------------------------------------------
% * compute
% -------------------------------------------

[ypdf,xpdf] = hist(vec,50);

xcdf = sort(vec);
ycdf = (1:length(xcdf)) ./ length(xcdf);

u = min(xcdf(ycdf >= (1-crit_p)));


% -------------------------------------------
% * plot
% -------------------------------------------
if doplot
    for i = 1:doplot
        if length(fh) >= i, figure(fh(i)), else figure, end
        hold on
        if i == 1, h(1) = plot(xpdf,ypdf./max(ypdf),'r','LineWidth',2);, end
        if i == 2, h(2) = plot(xcdf,ycdf,'r','LineWidth',2);, end
    end
end

return