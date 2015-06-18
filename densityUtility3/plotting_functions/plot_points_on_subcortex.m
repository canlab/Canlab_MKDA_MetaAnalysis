function [cl,han,surfhan] = plot_points_on_subcortex(DB,name,colors,condf, varargin)
% [cl,han,surfhan] = plot_points_on_subcortex(DB,name,colors,condf, varargin)
% 
% Extract and plot points for specific structure/structures.
%
% name input specifies structure to extract and plot on. Can be:
% - 'brainstem-thalamus', 'limbic', brainstem, thalamus, hypothalamus, caudate, putamen, or other names
%       recognized by addbrain
% - a mask image name of your choosing
% - already-defined clusters
%
% classes is
% a vector of integers (condf) or empty
%
% colors is
% a cell array with one color specification
% or one for each integer in classes
%
% condf is
% a vector of integers with points to plot in different colors; one
% integer per color entry
%
% textcodes is an optional input after condf
% a cell array of unique text labels to use instead of points
%
% an optional input after textcodes is an integer for how close is close_enough to target structure, in mm
%
% Examples:
% [ind,nms,condf] = string2indicator(DB.Mode);
% wh = [1 3];
% condf = indic2condf(ind(:,wh)); colors = {'ro' 'b^'};
% tor_fig; [cl,han] = plot_points_on_subcortex(DB,name,colors,condf)
%
% DB has xyz field with coordinates
% DB.Contrast is a vector of integers
% DB.textcodes is a cell array of text labels
% [cl, pointhan, surfhan] = plot_points_on_subcortex(DB,'limbic', {'ro'}, DB.Contrast, DB.textcodes);
% [cl, pointhan, surfhan] = plot_points_on_subcortex(DB,'limbic', {'ro'}, DB.Contrast);

global textcodes  % i used a global because i'm lazy...
global close_enough

textcodes = [];
close_enough = 2;

if length(varargin) > 0
    textcodes = varargin{1};
end

if length(varargin) > 1
    close_enough = varargin{2};
end

hold on;

switch name
    case 'brainstem-thalamus'
        [cl,han{1},surfhan{1}] = show_cl(DB,'brainstem',colors,condf);
        [c,h,surfhan{end+1}] = show_cl(DB,'thalamus',colors,condf);
        cl = merge_clusters(cl,c);
        han{end+1} = h;
        
        [c,h,surfhan{end+1}] = show_cl(DB,'hypothalamus',colors,condf);
        cl = merge_clusters(cl,c);
        han{end+1} = h;
    case 'limbic'
        [cl,han{1},surfhan{1}] = show_cl(DB,'brainstem',colors,condf);
        [c,h,surfhan{end+1}] = show_cl(DB,'thalamus',colors,condf);
        cl = merge_clusters(cl,c);
        han{end+1} = h;
        
        [c,h,surfhan{end+1}] = show_cl(DB,'hypothalamus',colors,condf);
        cl = merge_clusters(cl,c);
        han{end+1} = h;
        
        [c,h,surfhan{end+1}] = show_cl(DB,'caudate',colors,condf);
        cl = merge_clusters(cl,c);
        han{end+1} = h;
        
        [c,h,surfhan{end+1}] = show_cl(DB,'amygdala',colors,condf);
        cl = merge_clusters(cl,c);
        han{end+1} = h;
        
        [c,h,surfhan{end+1}] = show_cl(DB,'nucleus accumbens',colors,condf);
        cl = merge_clusters(cl,c);
        han{end+1} = h;
    otherwise
        [cl,han,surfhan{1}] = show_cl(DB,name,colors,condf);
        
end

axis image; axis off; axis vis3d; view(90,10); lighting gouraud
scn_export_papersetup(400);

end % main function


function [cl,han,surfhan] = show_cl(DB,name,colors,condf)

global textcodes
global close_enough

surfhan = [];
han = [];

if isstr(name)
    switch name
        case 'brainstem'
            pname = which('spm2_brainstem.img');
            surfhan = addbrain(name);
        case 'thalamus'
            pname = which('spm2_thal.img');
            surfhan = addbrain(name);
        case 'hypothalamus'
            pname = which('spm2_hythal.img');
            surfhan = addbrain(name);
        case 'caudate'
            pname = which('spm2_caudate.img');
            surfhan = addbrain(name);
        case 'amygdala'
            load amy_clusters
            cl = amy;
            surfhan = addbrain(name);
    
        case {'hipp', 'hippocampus'}
            pname = which('spm2_hipp.img');
            surfhan = addbrain(name);
            
        case 'brainstem'

            pname = 'spm2_brainstem.mat';
            surfhan = addbrain(name);
        
        case 'pag'
            pname = which('spm5_pag.img');
            surfhan = addbrain(name);
            
        case 'left'
            pname = which('spm2_left.img');
            surfhan = addbrain(name);
            
        case 'right'
            pname = which('spm2_right.img');
            surfhan = addbrain(name);
            
        case 'nucleus accumbens'
            P = which('NucAccumb_clusters.mat');
            load(P)
            cl = cl(1:2);
            surfhan = addbrain(name);
            
            cl = rmfield(cl, 'descrip');  % so we can use descrip field later if data is entered there
            
        case 'putamen'
            P = which('carmack_more_clusters.mat'); load(P)
            cl = put;
            surfhan= addbrain(name);
        otherwise
            pname = name;
    end
    
    if ~exist('cl','var')
        cl = mask2clusters(pname);
    end
    
elseif isstruct(name)
    cl = name;
end

% get rid of small clusters
wh = cat(1,cl.numVox);
wh = find(wh < 50);
cl(wh) = [];

han = [];
if isempty(cl)
    disp('No clusters > 50 voxels');
else
    
    cl = database2clusters(DB, cl, close_enough);
    
    for i = 1:length(cl)
        
        if isempty(textcodes)
            mytextcodes = [];
        else
            mytextcodes = textcodes(cl(i).wh_points);
        end
        
        mycondf = condf(cl(i).wh_points);
        mycolors = colors(unique(mycondf));
        
        h = plot_points_on_brain(cl(i).XYZmm', mycolors, mycondf, mytextcodes);
        han = [han h];
    end
    
    drawnow
end

end % show_cl function

