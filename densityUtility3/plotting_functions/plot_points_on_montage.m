function [newax, pointhandles] = plot_points_on_montage(xyz, varargin)
% [newax, pointhandles] = plot_points_on_montage(xyz, varargin)
%
% Plots points on montage of slices
% - Solid brain slices or contour outlines
% - Points or text labels or both
% - Flexible slice spacing, colors, marker sizes/styles, axis layout (one row/standard square)
% - axial or saggital orientation
% - Multiple different sets of points can be plotted in different colors/text labels
% - if 'contrast' vector is entered, will average nearby points (12 mm) within same contrast
%
% Takes all inputs of plot_points_on_slice.  See help for additional
% documentation of options.  In addition:
% 'onerow' : arrange axes in one row
% 'slice_range' : [min max] values in mm for slices to plot
% 'colorcond' or 'condf' : vector of integers that defines colors to use
%                   if used, colors should be a cell array
% 'color' : cell array of colors for each unique (ordered) value of condf
% 'contrast' : unique contrast numbers for each set of points, integers  
%              - will average nearby points (12 mm) within same contrast
%
% also valid:
% 'marker', 'points', 'MarkerSize', etc.
% see help plot_points_on_slice.m
%
% Examples:
%
% plot_points_on_montage(DB.xyz)
% plot_points_on_montage(xyz, 'text', DB.textcodes);
% plot_points_on_montage(xyz, 'text', DB.textcodes, 'color', [.2 .2 1], 'onerow');
% plot_points_on_montage(xyz, 'text', DB.textcodes, 'color', [.2 .2 1], 'onerow', 'sagittal');
% plot_points_on_montage(xyz, 'text', DB.textcodes, 'color', [.2 .2 1], 'sagittal', 'slice_range', [-50 50]);
% plot_points_on_montage(xyz, 'text', DB.textcodes, 'color', [.2 .2 1], 'sagittal', 'slice_range', [-40 40], 'solid', 'spacing', 20, 'close_enough', 10, 'onerow');

pointhandles = {};

myview = 'axial';
doonerow = 0;
spacing = 8; % slice spacing, in mm
ovl = which('SPM8_colin27T1_seg.img');  % which('scalped_avg152T1.img');
textcodes = [];
texthandles = [];
slice_range = 'auto';
condf = [];
contrastvals = [];

color = 'k';
%     facecolor = 'k';
%     marker = 'o';
%     drawslice = 1;
%     markersize = 12;
%     orientation = 'axial';
%     disptype = 'contour';


    % ------------------------------------------------------
    % parse inputs
    % ------------------------------------------------------

    for i = 1:length(varargin)
        if ischar(varargin{i})
            switch varargin{i}
                % reserved keywords
                
                case {'onerow'}, doonerow = 1;
                 
                case 'slice_range', slice_range = varargin{i + 1};
                  
                case {'spacing'}, spacing = varargin{i+1}; 

                case {'text', 'textcodes'} % do not pass on to slice plot... 
                    textcodes = varargin{i + 1};
                    varargin{i+1} = [];
                    varargin{i} = [];
                     
                case {'sag', 'sagg','saggital','sagittal'}, myview = 'sagittal';
    
                case {'condf' 'colorcond'}, condf = varargin{i + 1};
                    
                % From plot_points_on_slice
                
%                 case {'noslice', 'nodraw'}, drawslice = 0;

                 case 'color' % do not pass on...
                     color = varargin{i+1}; 
                     varargin{i+1} = [];
                     varargin{i} = [];
                     
                case 'contrast'
                    contrastvals = varargin{i+1};
                    varargin{i+1} = [];
                    varargin{i} = [];
%                     
% These are redundant here, but will be passed into plot_points_on_slice
% and are OK.
%                 case 'marker', marker = varargin{i+1}; 
% 
%                 case {'slice', 'wh_slice'}, wh_slice = varargin{i+1};
% 
%                     
                 
% 
%                 case {'MarkerSize', 'markersize'}, markersize = varargin{i+1};
% 
%                 case {'MarkerFaceColor', 'markerfacecolor'}, facecolor = varargin{i+1};
% 
%                 case 'solid', disptype = 'solid';
% 
%                 case 'overlay', ovl = varargin{i + 1};
                                        
                %otherwise, warning(['Unknown input string option:' varargin{i}]);
            end
        end
    end
    
% SETUP
% -----------------------------------------------

setup_axes();

if isempty(condf)
    condf = ones(size(xyz, 1), 1);
    
    if iscell(color)
        error('Color should not be a cell array.');
    end
    color = {color};

else
    if ~iscell(color)
        error('When entering condf, enter colors cell with same number of entries.');
    end
end

if ~isempty(contrastvals)
    [contrastvals, xyz, color, condf, textcodes] = setup_contrasts_and_average(xyz, color, contrastvals, textcodes, condf, varargin);
end

u = unique(condf);
n = length(u);

% Do the work for each slice
% -----------------------------------------------

figure(slices_fig_h)

for i = 1:length(slice_vox_coords)
    axes(newax(i));
    
    drawstr = 'draw';
    
    % Dec 2013: Tor modified to handle conditions (condf) in plot_points_on_slice
    
    if ~isempty(textcodes)
        pointhandles{i} = plot_points_on_slice(xyz, 'color', color, 'MarkerFaceColor', color, myview, drawstr, 'text', textcodes, 'slice', slice_vox_coords(i), varargin{:});
        
        title(sprintf('%s%3.0f', textbase, slice_mm_coords(i)));
        
        delete(pointhandles{i});
        
    else
        pointhandles{i} = plot_points_on_slice(xyz, 'color', color, 'MarkerFaceColor', color, myview, drawstr, 'slice', slice_vox_coords(i), varargin{:});
    end
    
    drawstr = 'nodraw';
    
%     for j = 1:n % For each color code
%         
%         wh = condf == u(j);
%     
%         if ~isempty(textcodes)
%             pointhandles{i} = plot_points_on_slice(xyz(wh, :), 'color', color{j}, 'MarkerFaceColor', color{j}, myview, drawstr, 'text', textcodes(wh), 'slice', slice_vox_coords(i), varargin{:});
%             
%             title(sprintf('%s%3.0f', textbase, slice_mm_coords(i)));
%             
%             delete(pointhandles{i});
%             
%         else
%             pointhandles{i} = plot_points_on_slice(xyz(wh, :), 'color', color{j}, 'MarkerFaceColor', color{j}, myview, drawstr, 'slice', slice_vox_coords(i), varargin{:});
%         end
%         
%         drawstr = 'nodraw';
%         
%     end
    
end

% not necessary?
%equalize_axes(newax);




% ------------------------------------------------------
% INLINE FUNCTIONS
% ------------------------------------------------------

function setup_axes

overlay = which('SPM8_colin27T1_seg.img');
[volInfo, dat] = iimg_read_img(overlay, 2);

myviews = {'axial' 'coronal' 'sagittal'};   % for selecting SPM window
whview = find(strcmp(myviews, myview));
if isempty(whview), error('myview must be axial, coronal, or sagittal.'); end

switch whview
    case 1
        if strcmp(slice_range, 'auto')
            slice_range = [-45 70];
        end
        
        cen = [slice_range(1):spacing:slice_range(2)]';
        slice_mm_coords = cen;
        cen = [zeros(length(cen), 2) cen];
        xyzvox = mm2voxel(cen, volInfo.mat);
        slice_vox_coords = xyzvox(:, 3);
    case 2
        if strcmp(slice_range, 'auto')
            slice_range = [-100 65];
        end
        
        cen = [slice_range(1):spacing:slice_range(2)]';
        slice_mm_coords = cen;
        cen = [zeros(length(cen),1) cen zeros(length(cen),1)];
        xyzvox = mm2voxel(cen, volInfo.mat);
        slice_vox_coords = xyzvox(:, 2);
    case 3
        if strcmp(slice_range, 'auto')
            slice_range = [-70 70];
        end
        
        cen = [slice_range(1):spacing:slice_range(2)]';
        slice_mm_coords = cen;
        cen = [ cen zeros(length(cen), 2)];
        xyzvox = mm2voxel(cen, volInfo.mat);
        slice_vox_coords = xyzvox(:, 1);
end


myviews2 = {'sagittal' 'coronal' 'axial' };  % for selecting coord
whcoord = strmatch(myview, myviews2) ;

% get text string base
mystr = {'x = ' 'y = ' 'z = '};
textbase = mystr{whcoord};

% get optimal number of axes
num_axes = size(cen, 1);
rc = ceil(sqrt(num_axes));

slices_fig_h = figure; %create_figure(myview);
set(slices_fig_h, 'Color', 'w');

if doonerow
    ss = get(0, 'ScreenSize');
    set(gcf, 'Position', [round(ss(3)/12) round(ss(4)*.9) round(ss(3)*.9) round(ss(4)/7) ])
end

for i = 1:num_axes
    
    if doonerow
        newax(i) = subplot(1, num_axes, i);
    else
        newax(i) = subplot(rc, rc, i);
    end
    
    axis off;
end

end % setup axes

end % main function




%[contrastvals, xyz, color, condf, textcodes] = setup_contrasts_and_average(XYZ, color, condf, varargin);


function [contrastvals, xyz, color, condf, textcodes] = setup_contrasts_and_average(xyz, color, contrastvals, textcodes, condf, varargin)

% if length(varargin) < 3
%     disp(['Not plotting text codes.'])
%     textcodes = [];
%     
% elseif length(varargin) > 2
%     textcodes = varargin{3};
% end

% % make textcodes into cell, if not
% % This is for averaging nearby by textcode, if the study/contrast
% % grping is not a cell array of strings.
% if ~iscell(textcodes), tmp={};
%     for i = 1:size(textcodes,1), tmp{i} = num2str(textcodes(i));  end
%     textcodes = tmp';
% end
% 
% % if plotting text codes, make white
% if length(varargin) < 4 && ischar(color{1}) && ~isempty(textcodes)
%     for i = 1:length(color), color{i}(2) = '.'; end
% end

contrasttext = mat2cell(contrastvals, ones(size(contrastvals, 1), 1), 1);
for i = 1:length(contrastvals)
    contrasttext{i} = num2str(contrasttext{i});
end

disp(['Averaging nearby points within 12 mm with the same contrast value.']);

if size(xyz, 1) > 1 && ~isempty(contrastvals)
    
    
    % average nearby coordinates together!  12 mm
    n = size(xyz, 1);
    [xyz,contrasttext,order] = average_nearby_xyz(xyz,12,contrasttext);
    if ~isempty(contrastvals) && length(contrastvals) == n, contrastvals = contrastvals(order); end
    if ~isempty(textcodes) && length(textcodes) == n, textcodes = textcodes(order); end
    if ~isempty(condf) && length(condf) == n, condf = condf(order); end
    
    % second pass, 8 mm
    n = size(xyz, 1);
    [xyz,contrasttext,order] = average_nearby_xyz(xyz,8,contrasttext);
    if ~isempty(contrastvals) && length(contrastvals) == n, contrastvals = contrastvals(order); end
    if ~isempty(textcodes) && length(textcodes) == n, textcodes = textcodes(order); end
    if ~isempty(condf) && length(condf) == n, condf = condf(order); end

elseif ~isempty(textcodes)
    % one point; leave as-is
    
else
    textcodes = [];
    disp('Plotting all peaks; no averaging.');
end

end % setup contrasts...
