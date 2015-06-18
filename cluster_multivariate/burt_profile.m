function [v,x] = burt_profile(burt2,nms,n,varargin)
% [v,x] = burt_profile(burt2,nms,n)
% 
% v = values plotted
% x = scaled burt matrix
%
% given a matrix of cross-correspondences in two sets, with the first n
% variables being the first set, this function makes a line plot in which
% the correspondences between the first set (rows, line series) and the
% second set (columns, x-axis) are shown. 
%
% any input for 4th (final optional) argument plots on existing figure
%
% in the meta-analysis toolbox, set 1 should be task types, and set 2
% should be brain regions.
%
% SCALING of y-axis:
% this function scales the input burt2 so that all diagonals are one
% then, if the indicators that went into burt2 code 1 for a study having an activation, 
% the diagonals of the 1st set represent the number of studies of a
% particular type.  Dividing each row by this number, the plotted
% off-diagonal elements will this reflect the cross product of study (set1) x
% region (set2) divided by total number of studies, or the proportion of studies
% activating that region.
% This is reflected in the y-axis label.
%
% However, if the inputs you enter are not scaled in this way - e.g., rows
% are centered - then the y-axis scaling will be different.


if length(varargin) == 0
    figure('Color','w'), 
end

set(gca,'FontSize',18),hold on

x = diag(1./diag(burt2)) * burt2; % this normalizes cov by row var %corrcoef(imatx);

v = x(1:n,n+1:end)';
plot(v,'s-','LineWidth',2); 
legend(nms{1:n})
set(gca,'XTick',1:length(nms)-n)
set(gca,'XTickLabel',nms(n+1:end),'XLim',[.5 (length(nms)-n)+.5])
ylabel('% Studies activating') 
xlabel('Brain region')

drawnow