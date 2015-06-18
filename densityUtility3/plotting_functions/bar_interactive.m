function bar_interactive(images, xnames)
% bar_interactive(images, xnames)
%
% create interactive bar-plot that pops up in spm_orthviews window
%
% tor wager
%
% 2006.05.02 - Modified by Matthew Davidson


% find the spm window, or make one from a p-image
spm_handle = findobj('Tag','Graphics');
if isempty(spm_handle) || ~ishandle(spm_handle)
    cl = pmap_threshold;
    spm_handle = gcf();
end

% setup images
if ~exist('images', 'var') || isempty(images)
    images = spm_get(Inf,'*img','Select images for bar plot');
end
V_images = spm_vol(images);
data_images = spm_read_vols(V_images);


% setup tick names
if exist('xnames', 'var') && ~isempty(xnames)
    % we have it
    
elseif (~exist('xnames', 'var') || isempty(xnames)) && exist('Xinms.mat', 'file')
    load('Xinms');
    xnames = Xinms;
else
    xnames = [];
    %xnames = images;
    fprintf('xnames not defined. Not displaying tick names.\n');
end

% setup figure to draw to
barplot_handle = findobj('Tag', mfilename);
if(isempty(barplot_handle))
    barplot_handle = init_bar_window();
end

% set(gcf,'WindowButtonUpFcn','dat = bar_interactive_btnupfcn;')
set(spm_handle,'WindowButtonUpFcn', {@bar_interactive_callback, barplot_handle, V_images, data_images, xnames})

return


function barplot_handle = init_bar_window()
barplot_handle = tor_fig();
set(barplot_handle, 'Tag', mfilename);
%axes('Position',[.52 .08 .45 .38]);
return


% -------------------------------------------------------------------
% Callback
% -------------------------------------------------------------------

function bar_interactive_callback(spm_handle, event_data, barplot_handle, V_images, data_images, xnames) 
% the first two params (the source handle and event-related data) are
% required by all callbacks - see documentation
% they may be ignored if preferred; the event data may be empty



% if no images selected, don't display anything.
if isempty(V_images)
    fprintf('No images selected. Displaying nothing.\n');
    return
end

% activate window
if isempty(barplot_handle) || ~ishandle(barplot_handle)
    barplot_handle = init_bar_window();
else
    figure(barplot_handle);
end
set(gca,'FontSize',16);

% get coordinate
mm_coord = spm_orthviews('Pos');
vox_coord = mm2voxel(mm_coord',V_images(1));
mm_coord = round(mm_coord');

dat = squeeze(data_images(vox_coord(1),vox_coord(2),vox_coord(3),:));

hold off;
bar_handles = bar(dat);
set(bar_handles,'FaceColor',[.7 .7 .7]);

if ~isempty(xnames), set(gca,'XTickLabel',xnames);, end
title(sprintf('[x,y,z] = %3.0f, %3.0f, %3.0f',mm_coord(1),mm_coord(2),mm_coord(3)),'FontSize',16);


% table output
% ----------------------------------------
%Meta_interactive_table;
Meta_interactive_table_vox(0);

return