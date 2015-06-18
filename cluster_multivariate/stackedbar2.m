function stackedbar2(v,colors,labels)
% stackedbar(v,colors,labels)
% plots a stacked bar plot, so you can see proportions of each variable out
% of whole
% useful for plotting p(task|act)
% columns of v are tasks, rows are regions
% labels are labels for regions

figure('Color','w'); set(gca,'FontSize',18); hold on;

if isempty(colors), colors = {'b' 'g' 'r' 'y' 'c' 'm'};, end
v(isnan(v)) = 0;

h = barh(v,'stacked');

for i =1:length(h), set(h(i),'FaceColor',colors{i});,end

set(gca,'YTick',1:size(v,1),'YTickLabel',labels,'YLim',[.5 size(v,1) + .5])

return