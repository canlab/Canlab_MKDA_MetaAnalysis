function varargout = studyplot(Method,varargin)
% function varargout = studyplot(method,coordinate list file,argument for selection,plothandle (string))
%
% method 'slices' or 'slices2'
%	[dist,voxelcoord,plothandles,contourhandles] = studyplot('slices');
%
% method 'glass' or 'glassvert'
%	[dist,brainhandle,linehandles,sagghandle,axhandle] = studyplot('glass');
%   [dist] = studyplot('glass','workingmem_database.txt','cdef','Stimuli_gen');
%   [dist] = studyplot('glass','workingmem_database.txt','cdef','Stimuli_gen','select','strcmp(Paradigm,''sternberg'')');
%
% method 'surface', 'surface2', surface3'
%	surface is original surface (doesn't show all points), 2 is solid within glass, 3 is lateral views.
%	[dist,surfacecoords,distfromsurface,brainhandle,lighthandle,chosenvertices] = studyplot('surface');
%
% method 'medial': medial surfaces
% 
%
% inputs: 
%   first = Method, 2nd = coordfile, third - nth are optional arguments.
%           special: third argument = which slices (for slices2 plot)
%
% in any order:
%   'text'          use text string to plot instead of colored shapes
%   'cdef'          followed by variable in database to define colors by
%   'select'        followed by coordinate selection string
%   'handles'       followed by plot handles, to add to existing figures
%                   Handles should be row vector of at least as many handles as needed for plot.
% newregion     replace region values with this new cell vector of regions
%  'numbers'       plot study numbers instead of points
%
% Variable input arguments (always text label followed by input):
%    if strcmp(varargin{i-1},'text'), plottext = 1;,end
%    if strcmp(varargin{i-1},'contour'), contouron = 1;,end
%    if strcmp(varargin{i-1},'cblm'), special = 'cblm';,end
%    if strcmp(varargin{i-1},'bg'), special = 'bg';,end
%	if strcmp(varargin{i-1},'amy'), special = 'amy';,end
%    if strcmp(varargin{i-1},'dist'), dodistances = 1;,end
%    if strcmp(varargin{i-1},'newregion'), newregion = (varargin{i});,end
%    if strcmp(varargin{i-1},'numbers'), donumbers = 1;,else, donumbers = 0;,end
%   'cdef'     followed by variable in database to define colors by
%   'select'    followed by coordinate selection string
%
% Last edits: 8/20/02 / surface fig update 4/26/04 by Tor Wager
%
%examples:
%figure('Color','w'); studyplot('surface','emotrev8_emotion3_t.txt','~(strcmp(colors,''ko'')) & ~(strcmp(colors,''wo''))');
% figure('Color','w');
% studyplot('slices2',dbname,'numbers','cdef','Study');

surfacecutoff = 20;                 % plot points this distance or less from surface, in mm.
mymarkersize = 4;

if nargin > 1,
	coordfile = varargin{1};
	if isempty(coordfile),coordfile = 'coordlist.txt';,end
else coordfile = 'coordlist.txt';
end

plottext = 0; contouron = 0; dodistances = 0; plothandle = [];
special = [];

for i = 2:nargin
    if strcmp(varargin{i-1},'text'), plottext = 1;,end
    if strcmp(varargin{i-1},'contour'), contouron = 1;,end
    if strcmp(varargin{i-1},'cblm'), special = 'cblm';,end
    if strcmp(varargin{i-1},'bg'), special = 'bg';,end
	if strcmp(varargin{i-1},'amy'), special = 'amy';,end
    if strcmp(varargin{i-1},'dist'), dodistances = 1;,end
    if strcmp(varargin{i-1},'handles'), plothandle = (varargin{i});,end
    if strcmp(varargin{i-1},'newregion'), newregion = (varargin{i});,end
    if strcmp(varargin{i-1},'numbers'), donumbers = 1;,else, donumbers = 0;,end
    if strcmp(varargin{i-1},'cdef'), cdef = varargin{i};,end
    if strcmp(varargin{i-1},'select'), selectstring = varargin{i};,end
    
end

if ischar(coordfile)
    
% -------------------------------------------------------------
% * load coordinate file
% -------------------------------------------------------------
disp(['Loading coordinate file ' coordfile])
dbname = coordfile;
read_database

% OLD WAY OF READING DATABASE
%[study,scan,roi,coordsys,numsubj,gender,stimuli,method,emotion,valence,target,ref,other,cogtype,cogdemand,x,y,z,zscore,region,colors,mytext] ...
%    = textread(coordfile,'%s%s%s%s%n%s%s%s%s%s%s%s%s%s%s%n%n%n%n%s%s%s');
%[x,y,z,colors,study,task,imgmethod,zscore] = textread(coordfile,'%n%n%n%s%s%s%s%n');

if exist('newregion')==1, 
    region = newregion;
    disp('Region input variable found: replacing original regions in text file with new regions.')
    unique(region)
    if exist('ColumnNames') == 1
        
        if ~any(strcmp(ColumnNames,'region')), ColumnNames{end+1} = 'region';, end
    else
        warning('No ColumnNames.  Database not read correctly. Pausing here.')
        keyboard
    end
end

else 
    % -------------------------------------------------------------
    % * we have entered coordinates
    % -------------------------------------------------------------
    x = coordfile(:,1); y = coordfile(:,2); z = coordfile(:,3);
    
end
    
% -------------------------------------------------------------
% * Assign colors if none exist - new
% -------------------------------------------------------------
if ~(exist('colors') == 1) & exist('cdef')
    % mylevels = eval(['unique(' cdef ')']);
    if iscell(cdef)
        mylevels = unique(cdef);
    else
        mylevels = eval(['unique(' cdef ')']);
        cdef = eval(cdef);
    end
    disp(mylevels), whch = input('Enter vector with indexes of these levels to save in analysis (e.g., [1 2 4]): ');
    disp(['saving ' cell2mat(mylevels(whch)')])
    myc = input('Enter color definitions for each condition (e.g., rgby): ','s');
    myc = mat2cell(myc,1,ones(length(myc)));
    
    try
        colors = cell(length(x),1);
        colors(:) = num2cell(0);
    catch
        error('No x variable in database?  Coordinate vars must be named x,y,z.')
    end
    %myc = {'r' 'g' 'b' 'y' 'c' 'm' 'w'};
    for i = 1:length(whch)
        a = strcmp(cdef,mylevels{whch(i)});
        %eval(['a = strcmp(' cdef ',mylevels{whch(i)});'])
        colors(a) = mat2cell([myc{min(i,length(myc))} 'o']);
    end
    
    % -------------------------------------------------------------   
    % * eliminate all rows with empty colors from all variables!
    % ColumnNames is created by read_database
    % -------------------------------------------------------------
    a = str2mat(colors{:}); a = a(:,1);
    a = (a==0); 
    ColumnNames{end+1} = 'colors';
    
    for i = 1:length(ColumnNames)
        mycol = eval([ColumnNames{i}]);
        if length(mycol) == length(a)
            eval([ColumnNames{i} '(a) = [];'])
        else
            warning(['Variable ' ColumnNames{i} ' is wrong size! (length = ' num2str(length(mycol)) ')'])
        end
    end
    %disp(['Selected ' num2str(size(x,1)) ' points in selected levels of ' cdef])

elseif ~(exist('colors') == 1)
    error('No ''colors'' variable found; enter as column in database or define using ''cdef'' input argument.')
end


% -------------------------------------------------------------   
% * Select coordinates based on input selection argument
% -------------------------------------------------------------
if donumbers, mystudy = unique(study);, end



if exist('selectstring') == 1
    disp(['Selecting requested coordinates out of original ' num2str(size(x,1)) ' points.'])
    %try,eval(['x = x(' varargin{2} ',1);']),catch, varargin{2}, str=['x = x(' varargin{2} ',1);'],end
    %eval(['y = y(' varargin{2} ',1);'])
    %eval(['z = z(' varargin{2} ',1);'])
    %index = 1;
    %eval(['whichones = ' varargin{2} ';'])
    %for i = 1:size(colors,1)
    %    if whichones(i), 
    %        colors2{index,1} = colors{i,1};
    %        mytext2{index} = mytext{i};
    %        region2{index} = region{i};
    %        study2{index} = study{i};
    %        index = index + 1;,end
    %end
    %colors = colors2; mytext = mytext2; region = region2; study = study2;
    %clear colors2; clear mytext2; clear region2; clear study2;

    % new way
    
    eval(['a = ' selectstring ';'])
    a = ~a;     % reverse, so 1's should be eliminated
    for i = 1:length(ColumnNames)
        mycol = eval([ColumnNames{i}]);
        if length(mycol) == length(a)
            eval([ColumnNames{i} '(a) = [];'])
        else
            warning(['Variable ' ColumnNames{i} ' is wrong size! (length = ' num2str(length(mycol)) ')'])
        end
    end
    
end

% -------------------------------------------------------------   
% * Compute inter-point distances, if asked for
% -------------------------------------------------------------

if dodistances
disp(['Computing inter-point distances for ' num2str(size(x,1)) ' coordinates.']);drawnow
end
dist = [];

% figure distances between points
for i = 1:size(x,1)
	coords(i,:) = [x(i) y(i) z(i)];
	if dodistances
        firstc = [x(i) y(i) z(i)];
	    for j = 1:size(x,1)
		    secondc = [x(j) y(j) z(j)];
		    diff = firstc - secondc;
		    dist(i,j) = sqrt(diff(1)^2 + diff(2)^2 + diff(3)^2);
	    end
    end
end
varargout{1} = dist;

disp(['Selected ' num2str(size(coords,1)) ' points that meet conditions.'])




% -------------------------------------------------------------   
% * Plot everything
% -------------------------------------------------------------


switch Method
	
	
	
	
	
	
	
	
	
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%	
case 'slices',
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
disp('Starting slice routine.'); drawnow
% load in the canonical T1 and the origin
basename = 'single_subj_T1';
basename = 'scalped_single_subj_T1.img';
[array,hdr] = readim2(basename);
origin = hdr.origin(1:3,1)';

% array = permute(array,[3 2 1]);

% put in proper orientation - rotate image and origin
clear aa
for i = 1:size(array,3)
    aa(:,:,i) = rot90(array(:,:,i));
end
array = aa;
clear aa

% origin is spec from top, matlab plots from bottom.  reverse y 
origin(2) = size(array,1) - origin(2);               

% coordinates entered in mm should be scaled to voxels
coords(:,1) = round(coords(:,1) / hdr.xsize);
coords(:,2) = round(coords(:,2) / hdr.ysize);
coords(:,3) = round(coords(:,3) / hdr.zsize);

% reverse y for coordinates
coords(:,2) = -coords(:,2);
% express all coordinates in voxel space, converting from distance from the origin
for i = 1:size(coords,1)
coords(i,:) = coords(i,:) + origin;
end

% specify orientation and swap coordinates appropriately; rotate T1
switch orient
case 'sagg'
    array = rotim('ax2sagg',array);
    coords = rotim('ax2sagg',coords);
case 'cor'
    array = rotim('ax2cor',array);
    coords = rotim('ax2cor',coords);
end

% figure out range of slices to plot - lowest to highest
minz = min(coords(:,3));
maxz = max(coords(:,3));

% figure out how many windows in the figure and make the figure
      nslices = maxz - minz + 1;
      temp = factor(nslices);
      if size(temp,2) > 3
         ximgs = temp(1) * temp(2) * temp(3);
         yimgs = 1;
         for j = 4:size(temp,2)
            yimgs = yimgs * temp(j);
         end
	  elseif size(temp,2) > 2
  	 		ximgs = temp(1) * temp(2);
  	 		yimgs = temp(3);
	  elseif size(temp,2) == 2
  	 		ximgs = max(temp);
  	 		yimgs = min(temp);
      else 
         ximgs = 8;
         yimgs = 8;
     end

     
% image the slices
colormap(gray)
clim = [min(min(min(array))) max(max(max(array)))];
for i = minz:maxz
    h(i) = subplot(yimgs,ximgs,i-minz+1);               % store subplot handle so that element of h is the slice
    imagesc(array(:,:,i),[clim]);
    axis off
    axis image
    %j = get(gca,'Position'); j = [j(1:2) .16 .16];
    %set(gca,'Position', j);
end
drawnow
% plot the points with specified colors
for i = 1:size(coords,1)
    subplot(h(coords(i,3))); hold on;
    plot(coords(i,1),coords(i,2),colors{i},'MarkerSize',2)
    drawnow
end

% 3-D contour plot
figure;
hand = contourslice(array,[],[],coords(:,3)');
set(hand,'LineWidth',1.5)
set(gcf,'Color',[1 1 1])
view(3); colormap(gray); hold on; daspect([1 1 .1]); drawnow
for i = 1:size(coords,1)
    plot3(coords(i,1),coords(i,2),coords(i,3),colors{i},'LineWidth',4);
end

set(gca,'FontSize',16); xlabel('Left-Right'),ylabel('Post-Ant');zlabel('Inferior-Superior')
set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]); set(gca,'ZTickLabel',[])

varargout{2} = coords;
varargout{3} = h;
varargout{4} = hand;










%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%	
case 'glassvert',
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
disp('Initiating glass brain plot 1 of 1')
% 3-D volume plot
% PLOT 1 - glass brain, lines are projections
global faces; global vertices;
load surf_single_subj_T1

% --- make figure ----------------------------------------------
if ~isempty(plothandle), axes(plothandle(1))
else figure
end
hold on
% --- plot coords ----------------------------------------------
for i = 1:size(coords,1)
    ph(i,:) = plotvert(coords(i,:),colors{i},vertices,'y'); drawnow
end
hold on
% --- 3-D brain   ----------------------------------------------
p = patch('Faces',faces,'Vertices',vertices);
view(135,30)
set(p,'EdgeColor','none')
lighting gouraud; alpha(.3); lighth = lightangle(45,30);
set(p,'FaceColor',[.4 .4 .4])
daspect([1 1 1]);
set(gcf,'Color',[0 0 0])
axis on; grid on; drawnow


varargout{2} = p;
varargout{3} = ph;










%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%	
case 'glass',
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
global faces; global vertices;
load surf_single_subj_T1

if donumbers,
        disp('studies are, in numerical order:')
        mystudy
end
    
disp('Initiating glass brain plot 1 of 3')
% PLOT 2 - saggital glass brain
% --- make figure ----------------------------------------------
if ~isempty(plothandle), axes(plothandle(1))
    if size(plothandle,2) < 3,error(['Not enough plot handles: ' num2str(size(plothandle,2)) ' found.']),end
else figure
end
hold on
% --- plot coords ----------------------------------------------
for i = 1:size(coords,1)
    if donumbers,
        %disp('studies are, in numerical order:')
        %mystudy
        mytext = num2str(find(strcmp(mystudy,study{i})));
        disp(['Study ' study{i} ' number ' mytext ' coords ' num2str(coords(i,:)) ' region ' region{i}])
        text(coords(i,1),coords(i,2),coords(i,3),mytext,'Color',colors{i}(1),'FontWeight','b');
      else
          plot3(coords(i,1),coords(i,2),coords(i,3),colors{i},'MarkerSize',mymarkersize,'MarkerFaceColor',colors{i}(1));
      end
end

% --- 3-D brain   ----------------------------------------------
p2 = patch('Faces',faces,'Vertices',vertices,'EdgeColor','none','FaceColor',[.4 .4 .4]);
view(270,0)
lighting gouraud; alpha(.3); lighth = lightangle(45,30); grid on % axis off;
daspect([1 1 1]);
set(gcf,'Color',[1 1 1])
set(gca,'ZTick',-60:10:80)
set(gca,'ZLim',[-70 90])
set(gca,'YLim',[-110 80])
set(gca,'YTick',[-100:20:70])
set(gca,'FontSize',8)
drawnow

varargout{2} = p2;

disp('Initiating glass brain plot 2 of 3')
% PLOT 3 - axial glass brain
% --- make figure ----------------------------------------------
if ~isempty(plothandle), axes(plothandle(2))
else figure
end
hold on
% --- plot coords ----------------------------------------------
for i = 1:size(coords,1)
	if donumbers,
        mytext = num2str(find(strcmp(mystudy,study{i})));
        text(coords(i,1),coords(i,2),coords(i,3),mytext,'Color',colors{i}(1),'FontWeight','b');
      else
          plot3(coords(i,1),coords(i,2),coords(i,3),colors{i},'MarkerSize',mymarkersize,'MarkerFaceColor',colors{i}(1));
      end
end

% --- 3-D brain   ----------------------------------------------
p3 = patch('Faces',faces,'Vertices',vertices,'EdgeColor','none','FaceColor',[.4 .4 .4]);
view(0,90)
lighting gouraud; alpha(.3); lighth = lightangle(45,30); grid on% axis off;
daspect([1 1 1]);
set(gcf,'Color',[1 1 1])
set(gca,'YLim',[-110 80])
set(gca,'YTick',[-100:10:70])
set(gca,'XLim',[-80 80])
set(gca,'XTick',-70:20:70)
set(gca,'FontSize',8)
text(-70,60,80,'L','FontSize',24)
text(60,60,80,'R','FontSize',24)
camzoom(1.4)
drawnow

varargout{3} = p3;

disp('Initiating glass brain plot 3 of 3')
% PLOT 4 - coronal glass brain
% --- make figure ----------------------------------------------
if ~isempty(plothandle), axes(plothandle(3))
else figure
end
hold on
% --- plot coords ----------------------------------------------
for i = 1:size(coords,1)
	if donumbers,
        mytext = num2str(find(strcmp(mystudy,study{i})));
        text(coords(i,1),coords(i,2),coords(i,3),mytext,'Color',colors{i}(1),'FontWeight','b');
      else
          plot3(coords(i,1),coords(i,2),coords(i,3),colors{i},'MarkerSize',mymarkersize,'MarkerFaceColor',colors{i}(1));
      end
end

% --- 3-D brain   ----------------------------------------------
p4 = patch('Faces',faces,'Vertices',vertices,'EdgeColor','none','FaceColor',[.4 .4 .4]);
view(0,0)
lighting gouraud; alpha(.3); lighth = lightangle(45,30); grid on% axis off;
daspect([1 1 1]);
set(gca,'XLim',[-80 80])
set(gca,'ZTick',-60:10:80)
set(gca,'XTick',-70:20:70)
set(gca,'FontSize',8)
text(-70,60,80,'L','FontSize',24)
text(60,60,80,'R','FontSize',24)
camzoom(1.4)
camzoom(1.04)
drawnow

varargout{4} = p3;










%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%	
case 'surface',
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
xtractto = 1;		% num of mm to extract points to away from surface
cylinderrad = 49;	% square of radius to find points within for nearestVertex
svert = []; svertdist = [];

disp('Loading brain surface.')
global faces; global vertices;
%load surf_single_subj_T1
surffile = which('surf_brain_render_T1.mat'); load(surffile)

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

%plot3(0,0,0,'ys','MarkerSize',8,'MarkerFaceColor','y')
%for i = 0:15:360,view(i,30),drawnow,end

if exist('svert'),varargout{2} = svert;,end
if exist('svertdist'),varargout{3} = svertdist;,end
varargout{4} = p;
varargout{5} = lighth;
varargout{6} = outvert;

return

% OLD CODE - not used
disp('Initiating surface plot 2 of 2')
% PLOT surface of brain
figure; p = patch('Faces',faces,'Vertices',vertices);
view(135,30)
set(p,'EdgeColor','none')
lighting gouraud; alpha(1); lighth = lightangle(45,30); % axis off;
set(p,'FaceColor',[.7 .7 .7])
daspect([1 1 1]);
set(gcf,'Color',[1 1 1])
set(gca,'Position',[-.05 -.05 .775*1.5 .815*1.5])
hold on
drawnow

for i = 1:size(coords,1)
    if svertdist(i) <= 18,
		pt(i) = plot3(svert(i,1),svert(i,2),svert(i,3),colors{i},'LineWidth',2);
	else pt(i) = 0;,
	end
    linecol = colors{i}; linecol = linecol(1);
    ph(i,:) = plotvert(coords(i,:),colors{i},vertices,linecol);
    drawnow
end
varargout{2} = svert;
varargout{3} = svertdist;
varargout{4} = p;
varargout{5} = pi;












%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%	
case 'surface2',
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%

disp('Loading brain surface.')
global faces; global vertices;
%load surf_single_subj_T1
surffile = which('surf_brain_render_T1.mat'); load(surffile)

disp('Initiating surface plot 1 of 2')
% --- make figure ----------------------------------------------
if ~isempty(plothandle), axes(plothandle(1))
end
hold on


% plot points
plot3(0,0,0,'ys','MarkerSize',8,'MarkerFaceColor','y')
clear svert; index = 1;outvert = [0 0 0];
disp('  Preparing points to plot...')
svert = coords(i,:); svertdist = [];pt = [];
for i = 1:size(coords,1)
    
    % UPDATED 3/14/01 by Tor Wager - "brain within a brain"
    %[svert(i,:),svertdist(i)] = nearestVertex(coords(i,:),vertices);
    %fprintf('.')

	%if svertdist(i) <= surfacecutoff,
        % add 1 ( or - 1 ) to all dimensions and plot text.
    %        outvert(index,:) = [svert(i,1) svert(i,2) svert(i,3)]; index = index + 1;
    %        disp(['     Plotting point ' num2str(i) ' . Coord ' num2str(coords(i,1)) ' ' num2str(coords(i,2)) ' ' num2str(coords(i,3)) ' plotted as coord ' num2str(svert(i,1)) ' ' num2str(svert(i,2)) ' ' num2str(svert(i,3))])
            if plottext
                pt(i) = text(svert(i,1)+sign(svert(i,1)),svert(i,2)+sign(svert(i,2)),svert(i,3)+sign(svert(i,3)),...
                mytext{i},'Color',colors{i}(1),'FontSize',12,'FontWeight','bold');
            else    
                pt(i) = plot3(coords(i,1),coords(i,2),coords(i,3),colors{i},'MarkerSize',mymarkersize,'MarkerFaceColor',colors{i}(1));
                %pt(i) = plot3(svert(i,1)+sign(svert(i,1)),svert(i,2)+sign(svert(i,2)),svert(i,3)+sign(svert(i,3)),colors{i},'MarkerSize',mymarkersize,'MarkerFaceColor',colors{i}(1));
            end   
        %drawnow
        %else pt(i) = 0;,
        %end
end

%numtoplot = sum(svertdist <= surfacecutoff);
%disp(['  ...done! Chosen ' num2str(numtoplot) ' points within ' num2str(surfacecutoff) ' mm of surface.'])

% PLOT surface of outside brain
p = patch('Faces',faces,'Vertices',vertices);
view(135,30)
set(p,'EdgeColor','none')
lighting gouraud; alpha(.5); lighth = lightangle(45,30); % axis off;
camlight right;
set(p,'FaceColor',[.7 .7 .7])
daspect([1 1 1]);
set(gcf,'Color',[1 1 1]), grid on
%set(gca,'Position',[-.05 -.05 .775*1.5 .815*1.5])
axis off
camzoom(1.2)
drawnow

for i = 1:3 % x,y,z scaling factors for inner brain
    xdist = mean([abs(max(vertices(:,i))) abs(min(vertices(:,i)))]);
    xscale = (xdist - surfacecutoff) / xdist;
    vertices2(:,i) = vertices(:,i) * xscale;
end
p2 = patch('Faces',faces,'Vertices',vertices2,'EdgeColor','none','FaceAlpha',1,'FaceColor',[0 0 0]);



varargout{4} = svert;
varargout{3} = svertdist;
varargout{2} = [p p2];
varargout{5} = lighth;
varargout{6} = outvert;














%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%	
case 'surface3',
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%

disp('Loading brain surface.')
global faces; global vertices;
%load surf_single_subj_T1
surffile = which('surf_brain_render_T1.mat'); load(surffile)

disp('Initiating surface plot 1 of 2')
hold on

% ---------------------------------------------------------------
% * right lateral view *
% ---------------------------------------------------------------

% --- make figure ----------------------------------------------
if ~isempty(plothandle), axes(plothandle(1))
end
hold on

% ---- plot origin ---------------------------------------------
plot3(0,0,0,'ys','MarkerSize',8,'MarkerFaceColor','y')
clear svert; index = 1;outvert = [0 0 0];
disp('  Preparing points to plot...')

% ---- select based on laterality ------------------------------
pcoords = coords(coords(:,1)>0,:);		% x > 0
pcolors = colors(coords(:,1)>0,:);
whichcolors = unique(pcolors)			% list of unique colors

disp('	Finding distance from surface...')
% ---- select based distance to surface-------------------------
for i = 1:size(pcoords,1)
	[dummy,svertdist(i)] = nearestVertex(coords(i,:),vertices);
	if mod(i,50) == 0, fprintf('%d',i),end
end
pcoords = pcoords(svertdist <= surfacecutoff,:);		% x > 0
pcolors = pcolors(svertdist <= surfacecutoff,:);
		
% ---- add surfacecutoff to coords to make visible -------------
pcoords(:,1) = pcoords(:,1) + surfacecutoff + 100;

disp('	Plotting points...')
for i = 1:size(whichcolors,1)
	% ---- select based on color group -------------------------
	%plotcolors = pcolors(strcmp(pcolors,whichcolors{i}) == 1,:);
	svert = pcoords(strcmp(pcolors,whichcolors{i}) == 1,:);

	% ---- plot the points -------------------------------------
	hold on
	if plottext & ~isempty(svert)
		text(svert(:,1),svert(:,2),svert(:,3),...
                mytext{i},'Color',plotcolors{i},'FontSize',12,'FontWeight','bold');
	elseif ~isempty(svert)
		plot3(svert(:,1),svert(:,2),svert(:,3),whichcolors{i},'MarkerSize',mymarkersize,'MarkerFaceColor',whichcolors{i}(1));
	end
end
drawnow
disp(['  ...done! Chosen ' num2str(size(pcoords,1)) ' points to plot - but not all may be within ' num2str(surfacecutoff) ' mm and thus visible.'])

% PLOT surface of brain
% --------------------------------------------------------------
p = patch('Faces',faces,'Vertices',vertices);
view(90,0)
set(p,'EdgeColor','none')
lighting gouraud; alpha(1); lighth = lightangle(45,30); % axis off;
camlight right;
set(p,'FaceColor',[.7 .7 .7])
daspect([1 1 1]);
set(gcf,'Color',[1 1 1]), grid on
%set(gca,'Position',[-.05 -.05 .775*1.5 .815*1.5])
axis off
%camzoom(1.9)
camzoom(1.9)
camdolly(.3,0,0)
drawnow

% ---------------------------------------------------------------
% * left lateral view *
% ---------------------------------------------------------------

% --- make figure ----------------------------------------------
if ~isempty(plothandle), axes(plothandle(2))
end
hold on

% ---- plot origin ---------------------------------------------
plot3(0,0,0,'ys','MarkerSize',8,'MarkerFaceColor','y')
clear svert; index = 1;outvert = [0 0 0];
disp('  Preparing points to plot...')

% ---- select based on laterality ------------------------------
pcoords = coords(coords(:,1)<0,:);		% x < 0
pcolors = colors(coords(:,1)<0,:);
whichcolors = unique(pcolors)			% list of unique colors

disp('	Finding distance from surface...')
% ---- select based distance to surface-------------------------
for i = 1:size(pcoords,1)
	[dummy,svertdist(i)] = nearestVertex(coords(i,:),vertices);
	if mod(i,50) == 0, fprintf('%d',i),end
end
pcoords = pcoords(svertdist <= surfacecutoff,:);		% x > 0
pcolors = pcolors(svertdist <= surfacecutoff,:);
		
% ---- add surfacecutoff to coords to make visible -------------
pcoords(:,1) = pcoords(:,1) - surfacecutoff - 100;

disp('	Plotting points...')
for i = 1:size(whichcolors,1)
	% ---- select based on color group -------------------------
	plotcolors = pcolors(strcmp(pcolors,whichcolors{i}) == 1,:);
	svert = pcoords(strcmp(pcolors,whichcolors{i}) == 1,:);

	% ---- plot the points -------------------------------------
	hold on
	if plottext & ~isempty(svert)
		text(svert(:,1),svert(:,2),svert(:,3),...
                mytext{i},'Color',plotcolors{i},'FontSize',12,'FontWeight','bold');
	elseif ~isempty(svert)
		plot3(svert(:,1),svert(:,2),svert(:,3),whichcolors{i},'MarkerSize',mymarkersize,'MarkerFaceColor',whichcolors{i}(1));
	end
end
drawnow

disp(['  ...done! Chosen ' num2str(size(pcoords,1)) ' points to plot - but all may be within ' num2str(surfacecutoff) ' mm and thus visible.'])

% PLOT surface of brain
% --------------------------------------------------------------
p = patch('Faces',faces,'Vertices',vertices);
view(270,0)
set(p,'EdgeColor','none')
lighting gouraud; alpha(1); %lighth = lightangle(45,30); % axis off;
camlight(240,80)
lightangle(330,30)
lightangle(210,-20)
set(p,'FaceColor',[.7 .7 .7])
daspect([1 1 1]);
set(gcf,'Color',[1 1 1]), grid on
%set(gca,'Position',[-.05 -.05 .775*1.5 .815*1.5])
axis off
camzoom(1.59)
drawnow

varargout{2} = svert;
varargout{3} = p;
varargout{4} = lighth;
varargout{5} = outvert;















%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%	
case 'medial',
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
disp('Loading canonical T1 image.');drawnow
basename = 'single_subj_T1';
basename = 'scalped_single_subj_t1';
[array,hdr] = readim2(basename);

if ~exist('coords'), coords = [0 0 0];,end
if ~exist('colors'), colors = {'ko'};,end

%origin = hdr.origin(1:3,1)';
clear D
for i = 1:size(array,3)
    E(:,:,i) = rot90(array(:,:,i));
end
% left is right, right is left (radiological) - so reverse x here only, not on slices...- NO, THIS MAKES EM BACKWARDS
%coords(:,1) = -coords(:,1);
% reverse y for coordinates 
coords(:,2) = -coords(:,2);
coordscopy = coords;


% =================================================================== 
% * Left Medial *
% =================================================================== 
disp('Plotting left medial coordinates.');drawnow

coords = coordscopy;
origin = hdr.origin(1:3,1)';
D = E;
% adjust origin because you removed voxels
origin(1) = origin(1) - 48;
% origin is spec from top, matlab plots from bottom.  reverse y 
origin(2) = size(E,1) - origin(2);
whos array
origin
disp('  Adjusting coordiates from mm to brain voxels');drawnow
% coordinates entered in mm should be scaled to voxels

coords(:,1) = round(coords(:,1) / hdr.xsize);
coords(:,2) = round(coords(:,2) / hdr.ysize);
coords(:,3) = round(coords(:,3) / hdr.zsize);

% express all coordinates in voxel space, converting from distance from the origin
for i = 1:size(coords,1)
    coords(i,:) = coords(i,:) + origin;
end


% --- make figure ----------------------------------------------
if ~isempty(plothandle), axes(plothandle(1))
else figure
end
hold on

origin
% plot origin
plot3([origin(1)-5 origin(1)-5],[origin(2) origin(2)],[1 size(array,3)],'k')
plot3([origin(1)-5 origin(1)-5],[1 size(array,2)],[origin(3) origin(3)],'k')
hold on; drawnow

% ---- select based on laterality ------------------------------
pcoords = coords(coordscopy(:,1) >= -(surfacecutoff) & coordscopy(:,1) <= 0,:);			
pcolors = colors(coordscopy(:,1) >= -(surfacecutoff) & coordscopy(:,1) <= 0,:);		% svert is coordinates to plot
whichcolors = unique(pcolors)			% list of unique colors
disp(['plotting ' num2str(size(pcoords,1)) ' points.'])


for i = 1:size(whichcolors,1)
	% ---- select based on color group -------------------------
	plotcolors = pcolors(strcmp(pcolors,whichcolors{i}) == 1,:);
	svert = pcoords(strcmp(pcolors,whichcolors{i}) == 1,:);

	% ---- plot the points -------------------------------------
	hold on
	if plottext & ~isempty(svert)
		text(zeros(size(svert,1),1),svert(:,2),svert(:,3),...
                mytext{i},'Color',plotcolors{i},'FontSize',12,'FontWeight','bold');
	elseif ~isempty(svert)
		plot3(zeros(size(svert,1),1),svert(:,2),svert(:,3),whichcolors{i},'MarkerSize',mymarkersize,'MarkerFaceColor',whichcolors{i}(1));
	end
end
drawnow

% ---- plot medial surface ----
D(:,1:48,:) = [];
D(:,:,77:end) = [];
disp('  Plotting surface');drawnow
p1 = patch(isosurface(D, 50),'FaceColor',[1,.75,.65],...,...
    'EdgeColor','none');
p2 = patch(isocaps(D, 50),'FaceColor','interp',...
    'EdgeColor','none');
view(269,5); axis tight; axis image; axis off
colormap(gray(100)); camzoom(1.35);
camlight left; camlight; camlight left; lighting gouraud
isonormals(D,p1)
drawnow
hold on


% =================================================================== 
% * Right Medial *
% =================================================================== 
disp('Plotting right medial coordinates.');drawnow

disp('Right Medial')
coords = coordscopy;

origin = hdr.origin(1:3,1)';
D = E;
% adjust origin because you removed voxels
%origin(1) = origin(1) + -48;
% origin is spec from top, matlab plots from bottom.  reverse y 
origin(2) = size(E,1) - origin(2);   

% coordinates entered in mm should be scaled to voxels
coords(:,1) = round(coords(:,1) / hdr.xsize);
coords(:,2) = round(coords(:,2) / hdr.ysize);
coords(:,3) = round(coords(:,3) / hdr.zsize);

disp('  Adjusting coordiates from mm to brain voxels');drawnow
% adjust coords because you removed voxels
%    coords(:,1) = coords(:,1) + -48;
% express all coordinates in voxel space, converting from distance from the origin
for i = 1:size(coords,1)
    coords(i,:) = coords(i,:) + origin;
end


% --- make figure ----------------------------------------------
if ~isempty(plothandle), axes(plothandle(2))
else figure
end
hold on

origin
% plot origin
plot3([origin(1)+5 origin(1)+5],[origin(2) origin(2)],[1 size(array,3)],'k')
plot3([origin(1)+5 origin(1)+5],[1 size(array,2)],[origin(3) origin(3)],'k')
hold on; drawnow

% ---- select based on laterality ------------------------------
pcoords = coords(coordscopy(:,1) <= (surfacecutoff) & coordscopy(:,1) >= 0,:);			
pcolors = colors(coordscopy(:,1) <= (surfacecutoff) & coordscopy(:,1) >= 0,:);		
whichcolors = unique(pcolors)			% list of unique colors
disp(['plotting ' num2str(size(pcoords,1)) ' points.'])


for i = 1:size(whichcolors,1)
	% ---- select based on color group -------------------------
	plotcolors = pcolors(strcmp(pcolors,whichcolors{i}) == 1,:);
	svert = pcoords(strcmp(pcolors,whichcolors{i}) == 1,:);

	% ---- plot the points -------------------------------------
	hold on
	if plottext & ~isempty(svert)
		text(ones(size(svert,1),1).*(origin(1)+1),svert(:,2),svert(:,3),...
                mytext{i},'Color',plotcolors{i},'FontSize',12,'FontWeight','bold');
	elseif ~isempty(svert)
		plot3(ones(size(svert,1),1).*(origin(1)+1),svert(:,2),svert(:,3),whichcolors{i},'MarkerSize',mymarkersize,'MarkerFaceColor',whichcolors{i}(1));
	end
end
drawnow



% ---- plot medial surface ----
D = E;
D(:,43:end,:) = [];
D(:,:,77:end) = [];
disp('  Plotting surface');drawnow
p1 = patch(isosurface(D, 50),'FaceColor',[1,.75,.65], ...
    'EdgeColor','none');
p2 = patch(isocaps(D, 50),'FaceColor','interp', ...
    'EdgeColor','none');
%view(75,10); 
view(90,5)
axis tight; axis image; axis off
colormap(gray(100));camzoom(1.35);
camlight; camlight right; lighting gouraud
isonormals(D,p1)
drawnow














%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%	
case 'slices2',
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
% load in the canonical T1 and the origin
basename = 'scalped_single_subj_T1';
[array,hdr] = readim2(basename);
origin = hdr.origin(1:3,1)';

% array = permute(array,[3 2 1]);

% put in proper orientation - rotate image and origin
clear aa
for i = 1:size(array,3)
    aa(:,:,i) = rot90(array(:,:,i));
end
array = aa;
clear aa

% origin is spec from top, matlab plots from bottom.  reverse y 
origin(2) = size(array,1) - origin(2);               

% coordinates entered in mm should be scaled to voxels
coords(:,1) = round(coords(:,1) / hdr.xsize);
coords(:,2) = round(coords(:,2) / hdr.ysize);
coords(:,3) = round(coords(:,3) / hdr.zsize);

% reverse y for coordinates
coords(:,2) = -coords(:,2);
% express all coordinates in voxel space, converting from distance from the origin
for i = 1:size(coords,1)
coords(i,:) = coords(i,:) + origin;
end

% specify orientation and swap coordinates appropriately; rotate T1
switch orient
case 'sagg'
    array = rotim('ax2sagg',array);
    coords = rotim('ax2sagg',coords);
case 'cor'
    array = rotim('ax2cor',array);
    coords = rotim('ax2cor',coords);
end

% figure out range of slices to plot - lowest to highest
minz = min(coords(:,3));
maxz = max(coords(:,3));
if nargin > 3 & ~(isstr(varargin{3})), whichslices = varargin{3};,else whichslices = minz:maxz;,end 

if strcmp(special,'cblm')
    whichslices = 7:2:37;
elseif strcmp(special,'bg')
    whichslices = 29:2:48;
elseif strcmp(special,'amy')
	whichslices = 27:3:38;
end

whichslices
% plot the anatomical slices
[dummy3,dummy,h,dummy2,rows,cols] = readim2(array,'p','ax',[],whichslices);
colormap(gray)

if strcmp(special,'cblm')
    for i = 1:size(h,2)
        axes(h(i))
        axis image
        axis([15 80 55 100])
        camzoom(.8)
    end
    drawnow
end

noplotcount = 0;

% plot the points with specified colors
for i = 1:size(coords,1)
    temp = whichslices((whichslices-coords(i,3)).^2 == min((whichslices-coords(i,3)).^2));
    myslice(i) = temp(end);         % if halfway between slices, choose top one.
    mydistinmm(i) = hdr.zsize * (coords(i,3) - myslice(i)); 
    
    if abs(mydistinmm(i)) < 5,
        subplot(h(whichslices == myslice(i))); hold on;
        if plottext
                pt(i) = text(svert(i,1)+sign(svert(i,1)),svert(i,2)+sign(svert(i,2)),...
                mytext{i},'Color',colors{i}(1),'FontSize',10,'FontWeight','bold');
        else    
                pt(i) = plot(coords(i,1),coords(i,2),colors{i},'MarkerSize',3,'MarkerEdgeColor','None','MarkerFaceColor',colors{i}(1));
        end      
    else disp(['Distance is > 5 mm (' num2str(mydistinmm(i)) ') from any plotted slice for point ' num2str(coords(i,:))])
        noplotcount = noplotcount+1;    
    end    
end
set(gcf,'Color',[1 1 1])

disp(['distances are ' num2str(mydistinmm)])
disp([num2str(noplotcount) ' points not plotted due to out of range.'])

h=get(gcf,'Children')
for i=1:length(h),axes(h(i)),camzoom(1.3),end
for i=1:length(h),axes(h(i)),camzoom(1.3),end
for i=1:length(h),axes(h(i)),camzoom(1.3),end

if contouron
% 3-D contour plot
figure;
hold on
for i = whichslices
    hand = patch(isocaps(array(:,:,1:i),60),'FaceColor','interp','EdgeColor','none'); 
end
drawnow
colormap(gray)
set(gcf,'Color',[1 1 1])
view(3); colormap(gray); hold on; daspect([1 1 .1]);
drawnow
for i = 1:size(coords,1)
    plot3(coords(i,1),coords(i,2),coords(i,3),colors{i},'MarkerSize',3,'MarkerEdgeColor','None','MarkerFaceColor',colors{i}(1));
end
set(gca,'ZLim',[whichslices(1)-1 whichslices(end)+1])
axis off
drawnow
alpha(.8);
end

varargout{2} = coords;
varargout{3} = h;
varargout{4} = noplotcount;


end 	% end switch
return
