function [handles, wh_slice, my_z] = plot_points_on_slice(xyz, varargin)
% [handles, wh_slice, my_coords, texthandles] = plot_points_on_slice(xyz, varargin)
%
% Usage:
% handles = plot_points_on_slice(xyz, 'nodraw', 'color', [0 0 1], 'marker', 'o','close_enough',8);
%
%     Optional inputs:
%     {'noslice', 'nodraw'}, drawslice = 0;
%
%     'color', color = varargin{i+1}; varargin{i+1} = [];
%
%     'marker', marker = varargin{i+1}; varargin{i + 1} = [];
%
%     {'slice', 'wh_slice'}, wh_slice = varargin{i+1};
%
%     {'close', 'closeenough', 'close_enough'}, close_enough = varargin{i+1};
%
%     {'sagg','saggital','sagittal'}, orientation = 'sagittal';
%
%     {'MarkerSize', 'markersize'}, markersize = varargin{i+1};
%
%     {'MarkerFaceColor', 'markerfacecolor'}, facecolor = varargin{i+1};
%
%     'solid', disptype = 'solid';
%
%     'overlay', ovl = varargin{i + 1};
%
%     {'text', 'textcodes'}, textcodes = varargin{i + 1}; varargin{i + 1} = [];
%
%     {'condf' 'colorcond'}, condf = varargin{i + 1};
%
%     'points', plot flat points, which is faster. default is 'spheres'
%
% sagittal plot:
% figure; [handles, wh_slice, my_z] = plot_points_on_slice(xyz,'slice',45,'sagg','close',8,'markersize',10);
%
% Tor Wager, March 2007
% Documentation not complete; refer to code for input options, etc.
% Updated Jan 2011
% Updated Feb 2013 - tor - default is spheres rather than points
% Updated Dec 2013 - tor - handle condition function
%
% wh_slice is in voxels
% close enough is in mm
%
% Example:
% create_figure('Slice view'); [handles, wh_slice, my_z] = plot_points_on_slice(DB.xyz(strcmp(DB.Valence, 'neg'), :), 'slice', 23, 'color', [0 0 1], 'marker', 'o','close_enough',8);
% [handles2, wh_slice, my_z] = plot_points_on_slice(DB.xyz(strcmp(DB.Valence, 'pos'), :), 'slice', 23, 'color', [1 1 0], 'marker', 'o','close_enough',8);
%
% Parent function example:
% newax2 = plot_points_on_montage(DB.xyz, 'condf', DB.condf, 'color', DB.color2, 'contrast', DB.contrast, 'axial', 'slice_range', [-40 40], 'solid', 'spacing', 10, 'close_enough', 8, 'onerow');
%
% More examples:
% plot_points_on_slice(PLOTINFO{1}.xyz, 'solid', 'wh_slice', 50, 'condf', PLOTINFO{1}.condf, 'MarkerSize', 6, 'color', PLOTINFO{1}.colors);
% plot_points_on_slice(PLOTINFO{1}.xyz, 'solid', 'wh_slice', 50, 'condf', PLOTINFO{1}.condf, 'MarkerSize', 6, 'color', PLOTINFO{1}.colors, 'points', 'MarkerFaceColor', PLOTINFO{1}.colors);
% plot_points_on_slice(PLOTINFO{1}.xyz, 'solid', 'wh_slice', 40, 'condf', PLOTINFO{1}.condf, 'MarkerSize', 6, 'color', PLOTINFO{1}.colors, 'MarkerFaceColor', PLOTINFO{1}.colors, 'saggital');

% defaults
% ------------------------------------------------------

color = 'k';
facecolor = 'k';
marker = 'o';
drawslice = 1;
close_enough = 10;   % in mm
markersize = 12;
orientation = 'axial';
disptype = 'contour';
ovl = which('SPM8_colin27T1_seg.img');  % which('scalped_avg152T1.img');
textcodes = [];
texthandles = [];
condf = []; % color codes
pointtype = 'spheres';
docondf = 0;

% ------------------------------------------------------
% parse inputs
% ------------------------------------------------------

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % reserved keywords
            case {'noslice', 'nodraw'}, drawslice = 0;
                
            case 'color', color = varargin{i+1}; varargin{i+1} = [];
                
            case 'marker', marker = varargin{i+1}; varargin{i + 1} = [];
                
            case {'slice', 'wh_slice'}, wh_slice = varargin{i+1};
                
            case {'spacing', 'close', 'closeenough', 'close_enough'}, close_enough = varargin{i+1};
                
            case {'sagg','saggital','sagittal'}, orientation = 'sagittal';
                
            case {'MarkerSize', 'markersize'}, markersize = varargin{i+1};
                
            case {'MarkerFaceColor', 'markerfacecolor'}, facecolor = varargin{i+1};
                
            case 'solid', disptype = 'solid';
                
            case 'overlay', ovl = varargin{i + 1};
                
            case {'text', 'textcodes'}, textcodes = varargin{i + 1}; varargin{i + 1} = [];
                
            case {'condf' 'colorcond'}, condf = varargin{i + 1}; docondf = 1;
                
            case {'slice_range' 'onerow' 'draw'}
                % also not used here, but passed in by default from calling
                % function plot_points_on_montage.m
                % so no warning...
                
            case 'axial'
                % default; not used
                
            case 'points'
                % plot points; default is now spheres
                pointtype = 'points';
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

% Check and set a few more things

if ~docondf
    condf = ones(size(xyz, 1), 1); 
    color = {color}; 
    facecolor = {facecolor};
else
    if ~iscell(facecolor), facecolor = repmat({facecolor}, 1, length(color)); end
end

% textcodes: enforce cell of strings
if ~iscell(textcodes)
    textcodes = mat2cell(textcodes, ones(size(textcodes, 1), 1), 1);
end

for i = 1:length(textcodes)
    if ~ischar(textcodes{i})
        textcodes{i} = num2str(textcodes{i});
    end
end

% ------------------------------------------------------
% get other info
% ------------------------------------------------------


ovlV = spm_vol(ovl);

if ~exist('wh_slice', 'var') || isempty(wh_slice)
    xyzvox = mm2voxel(xyz, ovlV.mat, 1);
    modes = mode(xyzvox, 1);
    wh_slice = modes(3);
    fprintf(1,'No slice entered. Picking slice with most points: z = %3.0f\n', wh_slice);
end

fprintf(1,'Plotting points within %3.0f mm of slice\n', close_enough);

% ------------------------------------------------------
% draw brain slice
% ------------------------------------------------------


if drawslice
    
    v = spm_read_vols(ovlV);
    
    mm = voxel2mm([1 1 1; size(v)]', ovlV.mat);
    xcoords = linspace(mm(1,1), mm(1,2), size(v, 1));
    ycoords = linspace(mm(2,1), mm(2,2), size(v, 2));
    %zcoords = linspace(mm(3,1), mm(3,2), size(v, 3));
    
    switch orientation
        case 'axial'
            [X, Y] = meshgrid(xcoords, ycoords); Z = v(:,:,wh_slice)';
            %figure;
            switch disptype
                case 'contour'
                    contour(X, Y, Z);
                case 'solid'
                    hh = imagesc(xcoords, ycoords, Z);
                    
                otherwise
                    error('unknown display type. check code.');
            end
            
            colormap(gray);
            
        case 'sagittal'
            zcoords = linspace(mm(3,1), mm(3,2), size(v, 3));
            [X, Y] = meshgrid(ycoords, zcoords); Z = squeeze(v(wh_slice,:,:))';
            %figure;
            switch disptype
                case 'contour'
                    h = contour(X, Y, Z);
                    
                case 'solid'
                    hh = imagesc(ycoords, zcoords, Z);
                    
                otherwise
                    error('unknown display type. check code.');
            end
            
            colormap(gray);
            h = findobj(gca,'Type','hggroup');
            set(h,'LineWidth',1.5);
            set(gca,'XDir','Reverse')
            
            
        otherwise, error('Unknown orientation.');
            
    end
    
    axis image
    
    switch disptype
        case 'solid'
            % set transparent value for clear axes
            myAlphaData = double(abs(Z) > 0);
            
            % If we set alphadata to clear for BG and axis color to none, we get clear
            % axes
            set(hh, 'AlphaDataMapping', 'scaled', 'AlphaData', myAlphaData)
            set(gca, 'Color', 'none')
    end
    
end

% ------------------------------------------------------
% put points on slice
% ------------------------------------------------------

switch orientation
    case 'sagittal'
        % get z coord in mm of this slice
        tmp_xyzmm = voxel2mm([wh_slice 1 1]', ovlV.mat);
        my_z = tmp_xyzmm(1);
        
        wh = find(abs(xyz(:,1) - my_z) <= close_enough);
        hold on;
        
%         if strcmp(pointtype, 'points')
%             handles = plot(xyz(wh,2), xyz(wh,3), marker, 'Color', color, 'MarkerFaceColor',facecolor, 'MarkerSize', markersize);
%             
%         elseif strcmp(pointtype, 'spheres')
%             xyztoplot = [xyz(wh,2), xyz(wh,3) zeros(length(wh), 1)];
%             handles = cluster_image_sphere(xyztoplot, 'color', color, 'radius', round(markersize/3)+1);
%         end
%         
%         if ~isempty(textcodes)
%             texthandles = plottext(xyz(wh, 2), xyz(wh, 3), textcodes(wh),  color, orientation);
%         end

        
                % set coords
        xtoplot = xyz(wh,2);
        ytoplot = xyz(wh,3);
        condftoplot = condf(wh); % ok even if no condf (see above)
        
        if ~isempty(textcodes)
            textcodestoplot = textcodes(wh);
        else
            textcodestoplot = [];
        end
            
        [handles, texthandles] = plot_coords(xtoplot, ytoplot, condftoplot, color, facecolor, pointtype, marker, markersize, textcodestoplot, orientation);
        
        
        
    case 'axial'
        
        % get z coord in mm of this slice
        tmp_xyzmm = voxel2mm([1 1 wh_slice]', ovlV.mat);
        my_z = tmp_xyzmm(3);
        
        wh = find(abs(xyz(:,3) - my_z) <= close_enough);
        hold on;
        
        % set coords
        xtoplot = xyz(wh,1);
        ytoplot = xyz(wh,2);
        condftoplot = condf(wh); % ok even if no condf (see above)
        
        if ~isempty(textcodes)
            textcodestoplot = textcodes(wh);
        else
            textcodestoplot = [];
        end
            
        [handles, texthandles] = plot_coords(xtoplot, ytoplot, condftoplot, color, facecolor, pointtype, marker, markersize, textcodestoplot, orientation);
        
        
        
        
    otherwise, error('Unknown orientation.');
        
end

if strcmp(disptype, 'solid')
    set(gca, 'YDir', 'normal');
    white_colormap;
end

% set sphere lighting
if strcmp(pointtype, 'spheres')
    lighting gouraud
    switch orientation
        case 'axial'
            [az, el] = view; el = el - 45;
            lh = lightangle(az, el);
        case 'sagittal'
            lightangle(0, -135);
            lightangle(0, -135);
            lightangle(0, -135);
    end
    set(handles, 'FaceAlpha', .6,'SpecularColorReflectance', .7);
end

end % function


function texthandles = plottext(x, y, textcodes, color, orientation)

% adjust for font placement
switch orientation
    case 'sagittal'
        xshift = 6;
    case 'axial'
        xshift = -6;
        
end


texthandles = [];
for j = 1:length(textcodes)
    % text labels
    
    
    if ischar(color)
        
        texthandles(j) = text(x(j) + xshift, y(j), textcodes{j},'Color',color,'FontSize',12,'FontWeight','bold');
        
    else
        texthandles(j) = text(x(j)  + xshift, y(j), textcodes{j},'Color',color,'FontSize',12,'FontWeight','bold');
        
    end
    
end

end % plottext


function white_colormap

cm = colormap;
wh = find(all(cm < .01, 2));
wh = all(cm < .01, 2);
cm(wh, :) = repmat([1 1 1], sum(wh), 1);
colormap(cm)

end


% function [color, facecolor] = parse_condf(condf, color, facecolor)
%
% origcolor = color;
% [color, facecolor] = deal(cell(length(condf), 1));
%
% u = unique(condf);
%
% for i = 1:length(u)
%
%     wh = condf == i;
%     color(wh) = origcolor(i);
%
%     if ischar(origcolor{i})
%         myc = origcolor{i}(1); % color = 1st (hopefully)
%         facecolor(wh) = {myc};
%     else
%         facecolor(wh) = origcolor(i);
%     end
%
% end
%
% end % subfunction

function [handles, texthandles] = plot_coords(xtoplot, ytoplot, condftoplot, color, facecolor, pointtype, marker, markersize, textcodes, orientation)

% xtoplot, ytoplot are in reference to figure dims, and are diferent for
% different orientations.

[handles, texthandles] = deal([]);

u = unique(condftoplot);

for i = 1:length(u)
    
    % set color
        colortoplot = color{i};
        facectoplot = facecolor{i};
        
        wh = condftoplot == u(i);
    
    if strcmp(pointtype, 'points')
        
        handles = plot(xtoplot(wh), ytoplot(wh), marker, 'Color', colortoplot, 'MarkerFaceColor',facectoplot, 'MarkerSize', markersize);
        
    elseif strcmp(pointtype, 'spheres')
        xyztoplot = [xtoplot(wh), ytoplot(wh) zeros(sum(wh), 1)];
        handles = cluster_image_sphere(xyztoplot, 'color', colortoplot, 'radius', round(markersize/3)+1);
        
    end
    
    if ~isempty(textcodes)
        texthandles = plottext(xtoplot(wh), ytoplot(wh), textcodes(wh), colortoplot, orientation);
    end
    
end

end % subfunction


