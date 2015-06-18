function h = plotbarcounts(countmat,numplots,testregion,varargin)
% function h = plotbarcounts(countmat,numplots,testregion,varargin)
%	case 'critical', critical = varargin{i+1};,index = index+2;
%	case 'title', mytitle = varargin{i+1};,i = i+2;
%	case 'legend', mylegend = varargin{i+1};,i = i+2;
%	case 'ylabel', myylabel = varargin{i+1};,i = i+2;
%	case 'colors', colors = varargin{i+1};,i = i+2;
%	case 'plothandle', plothandle = varargin{i+1};,index = index+2;
%   case 'noshade', shadeon = varargin{i+1};
% Must have testregion declared as global!

colsperplot = size(countmat,2) ./ numplots;
colindex = 1;

% defaults
% ===================================================
colors = {'ys-';'mo-';'rd--';'bx-';'g--^'};

%testregion = {'amy';'hipp';'sl';'hy';'thal'; ...
	%	'ofc';'mpfc';'scc';'acc';'pcc';'ins';'tp'; ...
		%'lpfc';'psmc';'tc';'pc';'oc'; ...
		%'bg';'cb';'bstem'};

critical = 0;
arealabels = 0;
shadeon = 0;
myylabel = ['<- R lat.  Left-Right counts  L lat. ->'];
mylegend = {'Happiness','Fear','Anger','Sadness','Disgust'};
mytitle = ['Laterality of counts across areas for specific emotions'];

shadeband = [.5 5.5 12.5 17.5];
shadecolor = [.5 .2 .2];

areatext = {'Limbic','Paralimbic','Uni/heteromodal','Other'};

%global laptopfigposition
%figposition = laptopfigposition;    
%figposition = [187   265   994   587]; %desktop
figposition = [122.0000  610.0000  921.0000  369.0000];  %laptop canonical
% optional args
% ===================================================
for i = 1:2:nargin-3
	switch varargin{i}
	case 'critical', critical = varargin{i+1};
	case 'title', mytitle = varargin{i+1};
	case 'legend', mylegend = varargin{i+1};
	case 'ylabel', myylabel = varargin{i+1};
	case 'colors', colors = varargin{i+1};
	case 'plothandle', plothandle = varargin{i+1};
	case 'overallcounts', ocounts = varargin{i+1};
	case 'arealabels',arealabels = varargin{i+1};
    case 'noshade', shadeon = varargin{i+1};
	end
end

% add overall counts to legend, if applicable
% ===================================================
if exist('ocounts')
	for i = 1:size(mylegend,2)
		mylegend{i} = [mylegend{i} ' (' num2str(ocounts(i)) ')           '];
	end
end

% data range
% ===================================================
mymin = min(min(countmat)); mymin = mymin(1) - .08 * mymin(1);
mymax = max(max(countmat)); mymax = mymax(1) + .08 * mymax(1);
overallyrange = [mymin mymax];
ymin = overallyrange(1);
ymax = overallyrange(2);

% plots
% ===================================================
for i = 1:numplots
	%lineindex = 1;
	if exist('plothandle'),
		axes(plothandle(i));
	elseif numplots > 1,h(i) = subplot(numplots,1,i);
	end
	
	%for j = colindex:colindex+colsperplot-1
		hold on; 
		bar(countmat,'grouped');
		colormap(lines(colsperplot));
		%lineindex = lineindex+1;
	%end
	%colindex = colindex + colsperplot;
	set(gca,'XTickLabel',testregion)
	set(gca,'XTick',1:size(testregion,1))
	axis([0 size(testregion,1)+1 overallyrange])
	
	if critical
		fill([0 0 21 21],[11 critical critical 11],[.2 .2 .2])
		alpha(.3)
		fill([0 0 21 21],[-11 -critical -critical -11],[.2 .2 .2])
		alpha(.3)
	end
	
    if shadeon
	    fill([shadeband(1) shadeband(2) shadeband(2) shadeband(1)],[ymin ymin ymax ymax],shadecolor,'FaceAlpha',.2)
	    fill([shadeband(3) shadeband(4) shadeband(4) shadeband(3)],[ymin ymin ymax ymax],shadecolor,'FaceAlpha',.2)
	    %fill([shadeband(5) shadeband(6) shadeband(6) shadeband(5)],[ymin ymin ymax ymax],shadecolor,'FaceAlpha',.2)
    end
	grid on
	
	if arealabels
		for j = 1:size(shadeband,2)
			text(shadeband(j),-12,areatext{j},'FontSize',10)
		end
	end
	
end

    ylabel(myylabel,'FontSize',14)
    xlabel('Brain area','FontSize',18)
	legend(mylegend)
	if exist('plothandle'),
		title(mytitle,'FontSize',18);
	else
		subplot(numplots,1,1); title(mytitle,'FontSize',18)
	end
	set(gcf,'Color',[1 1 1],'Position',figposition)
	
return