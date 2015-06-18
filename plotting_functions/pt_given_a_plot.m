function [ind,colors] = pt_given_a_plot(eff,p,legstr,varnamecode,regionnames)
% [indices, colors] = pt_given_a_plot(studycount,totalstudycount,code,varnamecode,regionnames)
%
% tmp =
% pt_given_a_plot(clnew.COUNTS.effpt_given_a,clnew.COUNTS.Ppt_given_a,levels,'PTplot_',clnew.COUNTS(1).clusternames);
%
% given counts, does barplot and CHI2 analysis
% easy to use with output from dbcluster_contrast_table.m
%
% perc :    counts or percent studies to plot; row is region, column is
%           condition
% code :    cell string of legend for conditions

colors = {'bo' 'gs' 'r^' 'yv' 'cs' 'mo'};

yoffset = 20;   % start plot at this + 50
xoffset = -.4;  % start x at 0 - this

format compact


if iscell(eff)
    eff = cat(1,eff{:});
end

if iscell(p)
    p = cat(1,p{:});
end

p = p < .05;

yvals = 50;
for i = 2:size(eff,2)
    yvals = [yvals yvals(end) + 5];
end
yvals = yvals + yoffset;

% get start coordinates x and y for each box
[x,y] = find(p);
eff = eff(find(p));
ind = x;

x = x + xoffset;
cols = colors(y);
y = yvals(y);

if isempty(y),
    % no results
    y = 70;
end

set(gca,'YLim',[0 max(y)+10])
fill([0 size(p,1) size(p,1) 0],[min(y)-5 min(y)-5 max(y)+10 max(y)+10],[.7 .7 .7],'FaceAlpha',.5)
title('Increase in likelihood of task given activation')

for i = 1:length(eff)
    
    % x start, y start, color, text string
    h(i) = drawbox(x(i),y(i),cols{i}(1),sprintf('%3.2f',eff(i)),sign(eff(i)));

end





return



function h1 = drawbox(xoffset,yoffset,color,txt,posneg);
dur = .8;    % how long box is in x
height = 4; % how high box is in y

x = [0 1 1 0]; x = x * dur + xoffset;
y = [0 0 1 1]; y = y * height + yoffset;

if posneg > 0
    h1 = fill(x,y,color,'FaceAlpha',.5);
else
    h1 = fill(x,y,[1 1 1]);
    set(h1,'EdgeColor',color,'LineWidth',2);
end

text(xoffset+.1*dur,yoffset+.5*height,txt,'FontSize',16,'FontWeight','b');

return





mytitle = [];
myylabel = ['Percentage of studies'];
mylegend = legstr;         %{'Spatial','Object','Verbal'};

for i = 1:length(mylegend), mylegend{i}(mylegend{i} == '_') = ' ';,end

disp('--- --- --- --- --- --- --- --- ---')

% make color map for bar graphs
% ======================================================
figure;
for i = 1:length(colors)
    h = plot(0,0,colors{i});
    cmap(i,:) =  get(h,'Color');
end
cmap = cmap(1:length(code),:);
close


fig2 = figure('Position',[122.0000  610.0000  921.0000  369.0000]);
set(gca,'FontSize',18)
myplot2 = subplot(1,1,1);
plotbarcounts(perc,1,regionnames,'title',mytitle,'ylabel',myylabel,'legend',mylegend,'plothandle',myplot2,'overallcounts',totalstudycount,'arealabels',0);
xlabel('')
%ylabel('% of total studies','FontSize',16)
colormap(cmap)
drawnow





