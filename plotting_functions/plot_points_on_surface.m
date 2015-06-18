function [svert,svertdist] = plot_points_on_surface(coords,varargin)
% function [svert,svertdist] = plot_points_on_surface(coords,colors)

plothandle = [];
surfacecutoff = 20;
plottext = 0;
mymarkersize = 8;

if length(varargin)>0, colors = varargin{1};, 
else
    for i = 1:size(coords,1)
        colors{i,1} = 'ro';
    end
end

if length(varargin)>1,
    plottext = 1;
    mytext = varargin{2};
end

xtractto = 1;		% num of mm to extract points to away from surface
cylinderrad = 49;	% square of radius to find points within for nearestVertex
svert = []; svertdist = [];

disp('Loading brain surface.')
global faces; global vertices;
%load surf_single_subj_T1
surffile = which('surf_brain_render_T1.mat'); 
surffile = which('surf_single_subj_T1_gray.mat');
%surffile = which('surf_Carmack_t1.mat');
load(surffile)

disp('Initiating surface plot 1 of 2')
% --- make figure ----------------------------------------------
if ~isempty(plothandle), axes(plothandle(1))
else figure
end
hold on


% plot points
plot3(0,0,0,'ys','MarkerSize',8,'MarkerFaceColor','y')
clear svert; index = 1;outvert = [0 0 0];
disp('  Preparing points to plot...')
for i = 1:size(coords,1)
    [svert(i,:),svertdist(i)] = nearestVertex(coords(i,:),vertices,cylinderrad);
    %fprintf('.')

	if svertdist(i) <= surfacecutoff,
        % add 1 ( or - 1 ) to all dimensions and plot text.
            outvert(index,:) = [svert(i,1) svert(i,2) svert(i,3)]; index = index + 1;
            disp(['     Plotting point ' num2str(i) ' . Coord ' num2str(coords(i,1)) ' ' num2str(coords(i,2)) ' ' num2str(coords(i,3)) ' plotted as coord ' num2str(svert(i,1)) ' ' num2str(svert(i,2)) ' ' num2str(svert(i,3))])
            if plottext
                pt(i) = text(svert(i,1)+xtractto*sign(svert(i,1)),svert(i,2)+xtractto*sign(svert(i,2)),svert(i,3)+xtractto*sign(svert(i,3)),...
                mytext{i},'Color',colors{i}(1),'FontSize',12,'FontWeight','bold');
            else    
                pt(i) = plot3(svert(i,1)+xtractto*sign(svert(i,1)),svert(i,2)+xtractto*sign(svert(i,2)),svert(i,3)+xtractto*sign(svert(i,3)),colors{i},'MarkerSize',mymarkersize,'MarkerFaceColor',colors{i}(1));
            end   
        %drawnow
	else pt(i) = 0;,
	end
end

numtoplot = sum(svertdist <= surfacecutoff);
disp(['  ...done! Chosen ' num2str(numtoplot) ' points within ' num2str(surfacecutoff) ' mm of surface.'])

% PLOT surface of brain
p = patch('Faces',faces,'Vertices',vertices);
view(135,30)
set(p,'EdgeColor','none')
lighting gouraud; alpha(1); lighth = lightangle(45,30); % axis off;
camlight right;
set(p,'FaceColor',[.7 .7 .7])
daspect([1 1 1]);
set(gcf,'Color',[0 0 0])
%set(gca,'Position',[-.05 -.05 .775*1.5 .815*1.5])
axis off
camzoom(1.3)
drawnow

end