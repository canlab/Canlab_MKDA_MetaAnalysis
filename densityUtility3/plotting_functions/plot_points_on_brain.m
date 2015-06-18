function h = plot_points_on_brain(XYZ,varargin)
% function handles = plot_points_on_brain(XYZ,varargin)
%
% This function plots a 3-column vector of xyz points (coordinates from
% studies) or text labels for each point.
%
% - option to plot on a glass brain.  Four different views are created in 2 figures.
% - same options as in plot_points_on_surface2.m, but does not extract
% (move) coordinates to surface
% - can also plot on existing surfaces/figures; used in
% plot_points_on_subcortex
%
% Usage:
% --------------------------------
% plot_points_on_brain(XYZ,[{color(s)}, colorclasses, {textmarkers/contrastcodes}, suppresstextdisplayflag, addbrainsurf])
%
% An optional 2nd argument is a color, or list of colors.
% An optional 3rd argument is a vector of integers to classify the points
% into groups.  Each group gets a color in the colors vector.
% An optional 4th argument is a cell array of text markers to plot
% instead of points. Points will be averaged within 12 mm if they have
% the same text label -- the idea is that text labels code unique
% contrasts within studies.
% 5th arg: do the averaging, but suppress text display; plot points
% instead.
%
% Note: You may not always be able to see text on the figure automatically;
% try addbrain
%
% examples:
% --------------------------------
% plot_points_on_brain(saddecrease,{'go'});
% plot_points_on_brain(sadpts,{'bo' 'rs' 'gv'},methi);
% plot_points_on_brain(saddecrease,{cell vec of all colors});
% plot_points_on_brain(saddecrease,{'go'},[],1);
%
% h = plot_points_on_brain(DB.xyz, {'bo'}, DB.Contrast, DB.textcodes);
%
% by tor wager
%
% Modified aug 06 by tor wager to add option to not add brain surface

% Modified Jan 11 by tor to add text label plotting.

h = [];
addbrainsurface = 0;


if size(XYZ, 2) ~= 3, error('XYZ must have n rows and 3 columns.'); end
if  isempty(XYZ)
    disp('No points to plot. Doing nothing.');
    return
end

% --------------------------------
% set up colors
% --------------------------------

[colors, classes] = setup_colors(XYZ, varargin{:});

% --------------------------------
% set up text
% --------------------------------

[textcodes, XYZ, colors, classes] = setup_text_and_average(XYZ, colors, classes, varargin{:});


%
%
%
%     colors = {'ro' 'go' 'bo' 'yo' 'co' 'mo' 'ko' 'r^' 'g^' 'b^' 'y^' 'c^' 'm^' 'k^'};
%
%     if length(varargin) > 0, colors = varargin{1}; end
%     if length(varargin) > 1, classes = varargin{2}; else classes = ones(size(XYZ,1),1); end
%     if length(varargin) > 2, addbrainsurface = varargin{3};  end
%
%     if isempty(classes), classes = ones(size(XYZ,1),1); end
%     if ~iscell(colors), colors = {colors}; end
%     if isempty(colors), colors = {'ro' 'go' 'bo' 'yo' 'co' 'mo' 'ko' 'r^' 'g^' 'b^' 'y^' 'c^' 'm^' 'k^'};, end
%
%     % if we enter a vector of colors for each point, set up classes and colors
%     if length(colors) == size(XYZ,1)
%         classes = zeros(size(colors));
%         coltmp = unique(colors);
%
%         for i = 1:length(coltmp)
%             wh = find(strcmp(colors,coltmp{i}));
%             classes(wh) = i;
%         end
%         colors = coltmp';
%     end

if addbrainsurface
    f1 = figure('Color','w'); set(gca,'FontSize',18),hold on
end

u = unique(classes);
if length(colors) < length(u)
    colors = repmat(colors, 1, length(u));
end


% Plot points and/or text
% ---------------------------------------------------------------
if isempty(textcodes)
    % No text codes - points only
    
    fprintf('Plotting %3.0f points\n', size(XYZ, 1));
    
    for clas = 1:length(u)
        
        whpoints = classes==u(clas);
        myxyz = [XYZ(whpoints, 1), XYZ(whpoints, 2), XYZ(whpoints, 3)];
        
        mycolor = colors{clas};
        mylinecolor = 'wo';
        myfacecolor = mycolor;
        
        if ischar(mycolor)
            myfacecolor = myfacecolor(1);
            if length(mycolor) > 1
                mylinecolor(2) = mycolor(2);
            end
        end
            
        h = plot3(myxyz(:, 1), myxyz(:, 2), myxyz(:, 3), ...
            mylinecolor, 'MarkerFaceColor', myfacecolor, 'MarkerSize', 8);
        
        hold on
    end
    
else  % We have text labels
    
    % Plot text labels
    fprintf('Plotting %3.0f text labels\n', size(XYZ, 1));
    fprintf(1, 'Class\tx\ty\tz\tCode\t\n')
     
    for j = 1:length(textcodes)
        % text labels
        myclass = classes(j);
        wh = u == myclass;
        
        myxyz = [XYZ(j,1),XYZ(j,2),XYZ(j,3)];
        
        if ischar(colors{wh})  % colors is a text string
            
            pt(j) = text(myxyz(1), myxyz(2), myxyz(3),...
                textcodes{j},'Color',colors{wh}(1),'FontSize',12,'FontWeight','bold');
            
        else % colors is a 3-element vector
            
            pt(j) = text(myxyz(1), myxyz(2), myxyz(3), ...
                textcodes{j},'Color',colors{wh},'FontSize',12,'FontWeight','bold');
             
        end
        
        % Print table along with points

        fprintf(1, '%3.0f\t%3.0f\t%3.0f\t%3.0f\t%s\n', myclass, myxyz(1), myxyz(2), myxyz(3), textcodes{j});
        
    end
    
end



if addbrainsurface
    
    addbrain;
    drawnow;
    
    % Add white filling to prevent see-through
    
    y = [-105 70 70 -105]; z = [-60 -60 78 78]; x = [0 0 0 0];
    hold on; hh = fill3(x,y,z,[.9 .9 .9]); set(hh,'EdgeColor','none','FaceAlpha',1)
    
    
    view(90,0)
    [az,el] = view;
    h = lightangle(az,el); set(gca,'FontSize',18)
    
    
    [hh1,hh2,hh3,hl,a1,a2,a3] = make_figure_into_orthviews;
    
end

end % main function





% ---------------------------------------------------------------------
% ---------------------------------------------------------------------

% Sub-functions

% ---------------------------------------------------------------------
% ---------------------------------------------------------------------


function [colors, classes] = setup_colors(XYZ, varargin)

colors = {'ro' 'go' 'bo' 'yo' 'co' 'mo' 'ko' 'r^' 'g^' 'b^' 'y^' 'c^' 'm^' 'k^'};

if length(varargin) > 0, colors = varargin{1};  end
if length(varargin) > 1, classes = varargin{2};  else classes = [];  end

%if isempty(classes), classes = ones(size(XYZ,1),1); end
if ~iscell(colors), colors = {colors}; end
if isempty(colors), colors = {'ro' 'go' 'bo' 'yo' 'co' 'mo' 'ko' 'r^' 'g^' 'b^' 'y^' 'c^' 'm^' 'k^'};  end

% if we have entered a vector of colors for each point, set up classes and colors

if length(colors) == size(XYZ,1) && isempty(classes)
    classes = zeros(size(colors));
    coltmp = unique(colors);
    
    for i = 1:length(coltmp)
        wh = find(strcmp(colors,coltmp{i}));
        classes(wh) = i;
    end
    colors = coltmp';
    
else % Do some checks
    nclasses = length(unique(classes(classes~=0)));
    
    if length(colors) == 1
        disp('Using single color for all point classes because you input a single color');
        colors = repmat(colors, 1, nclasses);
        
    elseif length(colors) < nclasses
        
        disp('There seem to be too few colors in your input. Using standard colors')
        colors = scn_standard_colors(length(unique(classes(classes~=0))));
    end
    
end

end


function [textcodes, XYZ, colors, classes] = setup_text_and_average(XYZ, colors, classes, varargin)

if length(varargin) > 3
    disp(['Not plotting text codes.'])
    textcodes = [];
    
elseif length(varargin) > 2
    textcodes = varargin{3};
end

% make textcodes into cell, if not
% This is for averaging nearby by textcode, if the study/contrast
% grping is not a cell array of strings.
if ~iscell(textcodes), tmp={};
    for i = 1:size(textcodes,1), tmp{i} = num2str(textcodes(i));  end
    textcodes = tmp';
end

% if plotting text codes, make white
if length(varargin) < 4 && ~isempty(textcodes)
    for i = 1:length(colors)
        if ischar(colors{i})
            colors{i}(2) = '.';
        else
            % no need to do anything... colors{i} = [1 1 1];
        end
    end
end

disp(['Averaging nearby points within 12 mm with the same text code.']);

if size(XYZ, 1) > 1 && ~isempty(textcodes)
    
    % average nearby coordinates together!  12 mm
    [XYZ,textcodes,order] = average_nearby_xyz(XYZ,12,textcodes);
    classes = classes(order);
    
    % second pass, 8 mm
    [XYZ,textcodes,order] = average_nearby_xyz(XYZ,8,textcodes);
    classes = classes(order);

elseif ~isempty(textcodes)
    % one point; leave as-is
    
else
    textcodes = [];
    disp('Plotting all peaks; no averaging.');
end



end


