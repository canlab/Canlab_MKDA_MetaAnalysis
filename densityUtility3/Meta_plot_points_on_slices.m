function Meta_plot_points_on_slices(DB, MC_Setup)
 % Meta_plot_points_on_slices(DB, MC_Setup)
 %
 % This function plots points on multiple slices.  It can do it either on
 % solid slices or "outline" contours.  Right now, it's hard-coded for
 % contours, but the main function it runs, plot_points_on_slice.m, has
 % input options for either, so this can be easily changed.
 %
 % If you're using the densityUtility3 toolbox, you need two structures
 % created with those tools:
 % DB, created using Meta_Setup.m, stored in the file SETUP.mat in your analysis directory
 % MC_Info, created using Meta_Activation_FWE.m, stored in the file MC_Info.mat in your analysis directory
 %
 % Once you run an analysis, you'll have both of those available, and can
 % run this function.
 %
 % NOTE: This function will color-code points according to the first
 % contrast you set up in Meta_Activation_FWE('setup'...).
 % If you did not specify contrasts, you will have to modify your MC_Setup
 % structure to add them:
 % MC_Setup.connames = {'EmoType'}; % this is the column in DB you want to use to color-code points on the plot
% MC_Setup.Xinms = {'Cog_Motor' 'Emotion'}; % these are the names of the
% fields in the column you specified above to identify the points of different types
%
 % However, you can also run this as a stand-alone function, as shown by
 % the example below:
 %
 % Example: Running this as a stand-alone function
 % You would create "dummy" DB and MC_Setup fields with all the inputs this
 % program needs, and pass those in
 %
% xyz = [18 8 8; -34 0 -6; 12 14 14; 18 8 6; 28 24 -2; 6 -16 6];
% %vec = [1 1 2 2 2]; % 1 early pain placebo reductions, 2 intxn with nalox during early
% DB.xyz = xyz; img = which('avg152T1.nii');  DB.maskV = spm_vol(img);
% MC_Setup.Xinms = {'er' 'inx'};             % arbitrary labels for colors
% DB.conditionlist = {'er' 'er' 'inx' 'inx' 'inx' 'inx'}'     % list of labels for each point (should be entries in MC_Setup.Xinms)
% MC_Setup.connames = {'conditionlist'}; % name of the column with labels
% DB.radius_mm = 10;  % close enough
%
% Tor Wager, Sept 1, 2009

orientflag = 'sagg'; 

DB.x = DB.xyz(:, 1); DB.y = DB.xyz(:, 2); DB.z = DB.xyz(:, 3);

 img = which('avg152T1.nii');   % overlay
spm_image('init', img);

Meta_interactive_point_slice_plot(MC_Setup, DB);

% montage

colors = {'b' 'r' 'g' 'y' 'm' 'k'}; % these match Meta_interactive_point_slice_plot

[axh, whsl] = setup_montage(DB.xyz, DB.maskV, orientflag);

condition_list = DB.(MC_Setup.connames{1}); % condition list
for i = 1:length(whsl)

    drawflag = 'draw';
    axes(axh(i))
    pointhandles = plot_points_on_slice(DB.xyz(strcmp(condition_list, MC_Setup.Xinms{1}), :), 'slice', whsl(i), 'color', colors{1}, 'marker', 'o', ...
        'close_enough',DB.radius_mm, 'markersize', 4, orientflag);

    for j = 2:length(MC_Setup.Xinms)

        if ~isempty(pointhandles), drawflag = 'nodraw'; end
        pth = plot_points_on_slice(DB.xyz(strcmp(condition_list, MC_Setup.Xinms{j}), :), 'slice', whsl(i), 'color', colors{j}, 'marker', 'o', ...
            'close_enough',DB.radius_mm, 'markersize', 4, orientflag);

        pointhandles = [pointhandles pth];

    end

end


end

function [axh, whsl] = setup_montage(XYZmm, V, orientflag)
    
    % how many slices, which ones

    XYZ = mm2voxel(XYZmm, V, 1)';    % 1 no re-ordering, allows repeats      
    
    switch orientflag
        case 'axial'
            whsl = unique(XYZ(3, :));               % which slices to show
        case 'sagg'
            whsl = unique(XYZ(1, :));               % which slices to show
        otherwise
            error('Orientflag, hard-coded in this function, must be either ''axial'' or ''sagg''')
    end
    
    nsl = length(whsl) + 1;
    rc = ceil(sqrt(nsl));                          % rows and columns for plot
    h = [];

    f1 = create_figure('SCNlab_Montage', rc, rc);
    colormap gray;
    set(f1, 'Color', [1 1 1], 'MenuBar', 'none')
    
    axh = findobj(f1, 'Type', 'axes');
    axh = sort(axh);
    
    % make axes look nice (this still needs work)
    for i = 1:length(axh)
        p = get(axh(i), 'Position');
        
        set(axh(i), 'Position', [p(1) - p(1)*.05 p(2)+p(2)*.05 p(3) - p(3)*.05 p(4)+p(4)*.05]);
        
    end
    
    if length(axh) > length(whsl)
        for i = length(whsl) + 1 : length(axh)
            delete(axh(i));
        end
    end
end