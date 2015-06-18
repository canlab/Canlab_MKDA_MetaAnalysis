function dat = bar_interactive_btnupfcn;
% [dat] = bar_interactive_btnupfcn;
%
% function for loading and plotting hewma data from a voxel.  must have declared globals f,f2,VOL,EXPT
% see hewma_timeseries_plot, the shell function.

global han      % handle of gui
global P        % names of images

global xnames   % x-axis names

global VOL      % volume info for images
% global VDAT     % volume data

% get image names
if isempty(P), 
    try load SETUP class_avg_images; P = class_avg_images;, catch, end
    
    if exist('SETUP.mat') == 2 & isempty(P)
        P = spm_get(Inf,'*img','Select Class Avg Images for bar plot.');
        class_avg_images = P;
        save SETUP -append class_avg_images
    end
end

if isempty(P),
    P = spm_get(Inf,'*img','Select Images for bar plot or none to skip.');
end


% if no images selected, don't display anything.
if isempty(P), return, end
    
P = check_valid_imagename(P);


% get axis tick names
if isempty(xnames),
    try load SETUP Xinms; xnames = Xinms;, catch, end
end

try, VOL = spm_vol(P);,
catch, 
    disp('Cannot find images to extract in setup, or error reading images:');
    disp(P)
    disp('try >>global P  in base workspace, then use spm_get to choose image names.');
    return
    %P = spm_get(Inf,'*img','Select Images for bar plot or none to skip.');
end

if isempty(P), return, end

VOL(1).M = VOL(1).mat;    

% activate window
if isempty(han) || ~ishandle(han), 
    fh = findobj('Tag','Graphics');
    data = guidata(fh);
    uicontrol(fh,'Style','Frame','Position',[.45 .02 .52 .5],...
        'BackgroundColor',[.8 .8 .8]);  % colored frame
    
    han = axes('Position',[.52 .08 .45 .38]);, 

else, 
    axes(han);, 
    fh = findobj('Tag','Graphics');
    data = guidata(fh);
end
set(gca,'FontSize',16);

% get coordinate
coord = spm_orthviews('Pos');
coord = mm2voxel(coord',VOL(1));
mm = round(spm_orthviews('Pos')');


% extract data and plot
if ~isfield(data,'VDAT') || isempty(data.VDAT), 
    VDAT = spm_read_vols(VOL);
    data.VDAT = VDAT;
    guidata(fh,data);
end

dat = squeeze(data.VDAT(coord(1),coord(2),coord(3),:));

hold off;
barh = bar(dat); set(barh,'FaceColor',[.7 .7 .7]);

set(gca,'XTickLabel',xnames);
title(sprintf('[x,y,z] = %3.0f, %3.0f, %3.0f',mm(1),mm(2),mm(3)),'FontSize',16);


% table output
% ----------------------------------------
Meta_interactive_table_vox(0);        
        




return