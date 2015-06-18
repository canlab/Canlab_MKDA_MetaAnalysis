function [h, pt, p] = plot_points_on_surface2(XYZ,varargin)
% function [axishan, pointhan, surfhan] = plot_points_on_surface2(XYZ,[{color(s)}, colorclasses, {textmarkers/contrastcodes}, varargin)
%
% This function plots a 3-column vector of xyz points (coordinates from
% studies) on a glass brain.  Four different views are created in 2
% figures.
%
% Optional inputs:
% ------------------------------------------------------------------------
% color: Cell array of colors, one per colorclass, or single color if all points should be same color
% colorclasses: a vector of integers to classify the points into color categories. Each group gets a color in the colors vector.
% textmarkers: a cell array of letter/number codes to use as markers
%  - peaks within 12 mm with the same text code are averaged together for clarity of display
%  - the idea is that textcodes identify unique contrasts from the same  study
%  - 'notext': suppresses the use of text labels, but leaves the averaging
%
% Notes:
% This function may not plot as many points as you enter, if they are a)
% not near enough the surface, b) contain exactly duplicated coordinates,
% or c) contain peaks within 12 mm with the same text code (treated as
% same-contrast).
% 
% Examples:
% ------------------------------------------------------------------------
% plot_points_on_surface2(XYZ,color);
% plot_points_on_surface2(XYZ,color,[],letter);    % plot letters!! (letter is cell array; averages nearby points with same letter code)
%
% Plot spheres, then set the material type on the spheres to 'default':
% [axishan, pointhan, surfhan] = plot_points_on_surface2(DB.xyz, DB.color, DB.condf);
% material(cat(2, pointhan{:}), 'default');
%
% Plot text codes:
% [axishan, pointhan, surfhan] = plot_points_on_surface2(DB.xyz, DB.color, DB.condf, DB.textcodes, 'nospheres');
%
% Suppress text codes and plot points:
% [axishan, pointhan, surfhan] = plot_points_on_surface2(DB.xyz, DB.color, DB.condf, DB.textcodes, 'notext', 'nospheres');
%
% Add text to existing plot, without re-rendering brain surfaces
% [axishan, pointhan, surfhan] = plot_points_on_surface2(DB.xyz, DB.color, DB.condf, DB.textcodes, 'nospheres', 'nobrains');
%
%
% Do nearby-averaging based on study ID; xyzem = [x y z]
% plot_points_on_surface2(xyzem,{'ro'},[],EMDB.Study,1);
%
% Select only some Emotions, and plot those in different colors
% [ind,nms,condf] = string2indicator(DB.Emotion);
% condf = indic2condf(ind(:,[2 3 4 5 7])); colors = {'ro' 'gv' 'md' 'ys' 'b^'};
% plot_points_on_surface2(DB.xyz,colors,condf,DB.Contrast,1);
% myp = findobj(gcf,'Type','Patch');
% set(myp,'FaceColor',[.85 .6 .5])
% for i = 1:6, subplot(3,2,i); axis off; end
%
% Make a legend for this figure
% nms = nms([2 3 4 5 7]);
% makelegend(nms,colors,1);
% scn_export_papersetup(200);
% saveas(gcf,'all_surf_emotion_pts_legend','png')
% 
% SEE ALSO: plot_points_on_subcortex.m, plot_points_on_brain.m

% Edited: Tor Wager, Jan 2011, for text plotting enhancement
% Dec 2012: Default behavior is spheres.  Use 'nospheres' to use points.
%
% Dec 2012: Tor : Fixed bug found by Luka with color assignment when some
% classes are empty in some plots. Affected spheres and points, but not
% text plots.  
%       - Fixed handle return

%f1 = figure('Color', 'w'); 


h = [];
pt = {};
p = [];

surfdist = 22;      % distance from surface to extract
dospheres = 1;
suppresstext = 0;
suppressbrains = 0;

if any(strcmp(varargin, 'nospheres')), dospheres = 0; end
if any(strcmp(varargin, 'suppresstext')), suppresstext = 1; end
if any(strcmp(varargin, 'notext')), suppresstext = 1; end
if any(strcmp(varargin, 'nobrains')), suppressbrains = 1; end

if suppressbrains
    % skip this
    f1 = gcf;
else
    f1 = create_figure('Surface Point Plot');
end

current_position = get(f1, 'Position');
set(f1, 'Position', [current_position(1:2) 1024 1280]);
set(gca,'FontSize',18);
hold on;


% --------------------------------
% set up colors
% --------------------------------

[colors, classes] = setup_colors(XYZ, varargin{:});

% --------------------------------
% set up text
% --------------------------------

[textcodes, XYZ, colors, classes] = setup_text_and_average(XYZ, colors, classes, varargin{:});

if suppresstext, textcodes = []; end
    
% --------------------------------
% make figures -- plot points
% --------------------------------

XYZ_orig = XYZ;     % original;
textcodes_orig = textcodes;
classes_orig = classes;

for i = 1:6 % for each view
    subplot(3,2,i);
    hold on
    
    h(i) = gca;
    
    XYZ = XYZ_orig;
    textcodes = textcodes_orig;
    classes = classes_orig;
    
    switch i
        case 1
            % top view, add to z
            wh = find(XYZ(:,3) <= 0);
            XYZ(wh,:) = [];
            XYZ(:,3) = XYZ(:,3) + surfdist;
            if ~isempty(textcodes), textcodes(wh) = []; end
            classes(wh) = [];
        case 2
            % bottom view
            wh = find(XYZ(:,3) >= 0);
            XYZ(wh,:) = [];
            XYZ(:,3) = XYZ(:,3) - surfdist;
            if ~isempty(textcodes), textcodes(wh) = []; end
            classes(wh) = [];
        case 3
            % left view
            wh = find(XYZ(:,1) >= 0);
            XYZ(wh,:) = [];
            XYZ(:,1) = XYZ(:,1) - surfdist;
            if ~isempty(textcodes), textcodes(wh) = []; end
            classes(wh) = [];
        case 4
             % right view
             wh = find(XYZ(:,1) <= 0);
            XYZ(wh,:) = [];
            XYZ(:,1) = XYZ(:,1) + surfdist;
            if ~isempty(textcodes), textcodes(wh) = []; end
            classes(wh) = [];
        case 5
            % front view
            %wh = find(XYZ(:,2) <= 0);
            %XYZ(wh,:) = [];
            %XYZ(:,2) = XYZ(:,2) + surfdist;
            %if ~isempty(textcodes), textcodes(wh) = []; end
            %classes(wh) = [];
            
            % right medial
            wh = find(~(XYZ(:,1) >= 0 & XYZ(:,1) <= 16));

             XYZ(wh,:) = [];
            XYZ(:,1) = XYZ(:,1) - 2.*surfdist;
            if ~isempty(textcodes), textcodes(wh) = []; end
            classes(wh) = [];
        case 6            
            % back view
            %wh = find(XYZ(:,2) >= 0);
            %XYZ(wh,:) = [];
            %XYZ(:,2) = XYZ(:,2) - surfdist;
            %if ~isempty(textcodes), textcodes(wh) = []; end
            %classes(wh) = [];
            
            % left medial
            wh = find(~(XYZ(:,1) <= 0 & XYZ(:,1) >= -16));
             XYZ(wh,:) = [];
            XYZ(:,1) = XYZ(:,1) + 2.*surfdist;
            if ~isempty(textcodes), textcodes(wh) = []; end
            classes(wh) = [];
    end
    
    hold on;
    if isempty(textcodes)
        % plot points, no text labels
        
        u = unique(classes);
        u(u == 0) = [];
        
        for clas = 1:length(u) % 1:max(classes)
            
            if ischar(colors{clas})
                % Text specification of color
                if dospheres
                    % SPHERES
                    pthan{clas} = cluster_image_sphere(XYZ(classes==u(clas), 1:3), 'color', colors{u(clas)}(1), 'radius', 4);
                    
                else
                    % POINTS
                    if length(colors{clas}) < 2
                        error('You must enter symbol and color, e.g., ''go''');
                    end
                    
                    pthan{clas} = plot3(XYZ(classes==u(clas),1),XYZ(classes==u(clas),2),XYZ(classes==u(clas),3), ...
                        ['w' colors{clas}(2)],'MarkerFaceColor',colors{u(clas)}(1),'MarkerSize',8);
                end
                
            else
                % 3-element color vector
                if dospheres
                    % SPHERES
                    pthan{clas} = cluster_image_sphere(XYZ(classes==u(clas), 1:3), 'color', colors{u(clas)}, 'radius', 4);
                    
                else
                    % POINTS
                    pthan{clas} = plot3(XYZ(classes==u(clas),1),XYZ(classes==u(clas),2),XYZ(classes==u(clas),3), ...
                        ['wo'],'MarkerFaceColor',colors{u(clas)},'MarkerSize',8);
                end
            end
            
        end % class
        
        pt{i} = cat(2, pthan{:});  % may break if empty?  fix if so...
        
    else
        % plot text labels instead
        
        u = unique(classes);
        u(u == 0) = [];
        
        for j = 1:length(textcodes)
            % text labels
            myclass = classes(j);
            wh = u == myclass;
            
            if ischar(colors{wh})
                
                pt{j} = text(XYZ(j,1),XYZ(j,2),XYZ(j,3),...
                    textcodes{j},'Color',colors{wh}(1),'FontSize',12,'FontWeight','bold');
                
            else
                pt{j} = text(XYZ(j,1),XYZ(j,2),XYZ(j,3),...
                    textcodes{j},'Color',colors{wh},'FontSize',12,'FontWeight','bold');
                
                
            end
            
        end
    end

    drawnow
end % for each view


% --------------------------------
% make figures -- add brains
% --------------------------------
if suppressbrains
    return
end

for i = 1:6 % for each view
    subplot(3,2,i);
    
    
    
        switch i
        case 1
            % top view, add to z
            p(i) = addbrain('hires'); set(p,'FaceAlpha',1); drawnow;
            view(0,90); [az,el] = view; 
            lh(i) = lightangle(az,el); 
            set(gca,'FontSize',18)
            material dull
        case 2
            % bottom view
            p(i) = addbrain('hires'); set(p,'FaceAlpha',1); drawnow;
            view(0,-90); [az,el] = view; 
%             lh(i) = lightangle(az,el); 
            set(gca,'FontSize',18)
            set(gca,'YDir','Reverse')
            material dull
            %lightRestoreSingle;
%             hh(1) = lightangle(90,0);
%             hh(2) = lightangle(0,0);
            
             hh = findobj(gca, 'Type', 'light');
             delete(hh)
            camlight(180, 0)


        case 3
            % left view
            p(i) = addbrain('hires'); set(p,'FaceAlpha',1); drawnow;
            view(270,0); [az,el] = view; 
            lh(i) = lightangle(az,el); 
            set(gca,'FontSize',18)
            material dull
        case 4
             % right view
             p(i) = addbrain('hires'); set(p,'FaceAlpha',1); drawnow;
            view(90,0); [az,el] = view; 
            lh(i) = lightangle(az,el); 
            set(gca,'FontSize',18)
            material dull
        case 5
            % front view
            %view(180,0); [az,el] = view; h(i) = lightangle(az,el); set(gca,'FontSize',18)

            % right medial
            p(i) = addbrain('hires right'); set(p,'FaceAlpha',1, 'FaceColor', [.5 .5 .5]); drawnow;
            view(270,0); [az,el] = view; 
            %lh(i) = lightangle(az,el); 
            set(gca,'FontSize',18)
            material dull
            set(p(i),'FaceColor',[.5 .5 .5])
            %camzoom(.9)
            lightRestoreSingle(gca);

        case 6            
            % back view
            %view(0,0); [az,el] = view; h(i) = lightangle(az,el); set(gca,'FontSize',18)
            
            % left medial
            view(90,0); [az,el] = view; 
            lh(i) = lightangle(az,el); set(gca,'FontSize',18)
            p(i) = addbrain('hires left'); 
            set(p(i),'FaceAlpha',1, 'FaceColor', [.5 .5 .5]); drawnow;
            material dull
            camzoom(1.04)
           
        end
    
        camzoom(1.3)
        axis image;
        axis off
        %lightRestoreSingle(gca);
        scn_export_papersetup(900);
        
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
if length(varargin) > 1, classes = varargin{2};  else classes = ones(size(XYZ,1),1);  end

if isempty(classes), classes = ones(size(XYZ,1),1); end
if ~iscell(colors), colors = {colors}; end
if isempty(colors), colors = {'ro' 'go' 'bo' 'yo' 'co' 'mo' 'ko' 'r^' 'g^' 'b^' 'y^' 'c^' 'm^' 'k^'};  end

% if we enter a vector of colors for each point, set up classes and colors

if length(colors) == size(XYZ,1)
    classes = zeros(size(colors));
    coltmp = unique(colors);
    
    for i = 1:length(coltmp)
        wh = find(strcmp(colors,coltmp{i}));
        classes(wh) = i;
    end
    colors = coltmp';
else
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

if length(varargin) < 3
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
if length(varargin) < 4 && ischar(colors{1}) && ~isempty(textcodes)
    for i = 1:length(colors), colors{i}(2) = '.'; end
end

if size(XYZ, 1) > 1 && ~isempty(textcodes)
    
    disp(['Averaging nearby points within 12 mm with the same text code.']);

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
