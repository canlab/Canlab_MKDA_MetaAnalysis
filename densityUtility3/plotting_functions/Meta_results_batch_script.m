%% CANlab MKDA meta-analysis results report
% "What's in the box": A series of tables and plots for a standard
% coordinate-based meta-analysis created with Meta_Activation_FWE in the
% CANlab MKDA meta-analysis toolbox.
%
% The script Meta_results_batch_script is meant to be run stand-alone or
% used in publish_meta_analysis_report to create an HTML report with
% results and plots.

% Display helper functions: Called by later scripts
% --------------------------------------------------------

dashes = '----------------------------------------------';
printstr = @(dashes) disp(dashes);
printhdr = @(str) fprintf('%s\n%s\n%s\n', dashes, str, dashes);

%% Results orthviews and tables with multi-threshold and subclusters
% Multi-threshold: Show results with voxel-wise FWER and 'stringent' 'medium' and
% 'lenient' cluster extent-corrected results in different colors.
%
% Tables show numbers for contiguous clusters and arrows for sub-peaks
% within each contiguous cluster using SPM software's local maximum-finding
% algorithm.

close all

cl = Meta_Activation_FWE('results', 1);

drawnow, snapnow

% Preparing regions so that higher thresholds are included in region maps
% for lower thresholds (this is not done in the code above).

h = fmri_data('Activation_FWE_height.img', 'noverbose');
e1 = fmri_data('Activation_FWE_extent_stringent.img', 'noverbose');
e2 = fmri_data('Activation_FWE_extent_medium.img', 'noverbose');
e3 = fmri_data('Activation_FWE_extent_lenient.img', 'noverbose');

mydat = h;
r_height = region(mydat);

mydat.dat = mydat.dat + e1.dat;
r_stringent = region(mydat);

mydat.dat = mydat.dat + e2.dat;
r_medium = region(mydat);

mydat.dat = mydat.dat + e3.dat;
r_lenient = region(mydat);

% Set up to extract
% Study proportion data for voxel-wise FWER-corrected regions
% First, we will load the MC_Info.mat file, which contains lots of
% information about the analysis.  The variable of main interest in this

% file is MC_Setup, which contains the contrast indicator maps.
% This is stored in MC_Setup.unweighted_study_data

load MC_Info

% MC_Setup =
%     unweighted_study_data: [231202×18 double]             % Voxels x contrast indicator maps (1/0)
%                   volInfo: [1×1 struct]                   % Info for mapping back into brain space
%                         n: [2 3 2 2 1 2 1 2 4 2 2 3 3 2 2 8 2 1] % Coordinates per study
%                       wts: [18×1 double]                  % Contrast weights
%                         r: 5                              % Radius in voxels
%                        cl: {1×18 cell}                    % contiguous blobs for each contrast, for Monte Carlo


%% Voxel-wise FWER: Table and slice montage with autolabeled regions

printhdr('Voxel-wise FWER');

if ~isempty(r_height)
    
    % warning off, r = cluster2region(cl{1}); warning on
    
    montage(r_height, 'colormap');
    create_figure('surface'); axis off; surface(r_height, 'cutaway');
    drawnow, snapnow
    
    [rpos, rneg] = table(r_height);
    montage([rpos rneg], 'regioncenters', 'colormap');
    
    drawnow, snapnow
    
    r_height = [rpos rneg];             % save with region names, also for data extraction below
    
    r_height = extract_data_and_plot(r_height, MC_Setup); % extract data from regions
    
else
    
    disp('No significant results.');
    
end


%% Stringent .001 cluster-extent corrected: Table and slice montage with autolabeled regions

printhdr('Voxel-wise FWER + stringent primary threshold (.001) cluster-extent');

if ~isempty(r_stringent)
    
    %r = cluster2region(cl{2});
    
    montage(r_stringent, 'colormap');
    create_figure('surface'); axis off; surface(r_stringent, 'cutaway');
    
    drawnow, snapnow
    
    [rpos, rneg] = table(r_stringent);
    montage([rpos rneg], 'regioncenters', 'colormap');
    
    drawnow, snapnow
    
    r_stringent = [rpos rneg];             % save with region names, also for data extraction below
    
    r_stringent = extract_data_and_plot(r_stringent, MC_Setup); % extract data from regions
    
else
    
    disp('No significant results.');
    
end


%% Medium .01 cluster-extent corrected: Table and slice montage with autolabeled regions
% Excluding voxels already reported at more stringent/localizable thresholds

printhdr('Voxel-wise FWER + medium primary threshold (.01) cluster-extent');

if ~isempty(r_medium)
    
    %warning off, r = cluster2region(cl{3}); warning on
    
    montage(r_medium, 'colormap');
    create_figure('surface'); axis off; surface(r_medium, 'cutaway');
    
    drawnow, snapnow
    
    [rpos, rneg] = table(r_medium);
    montage([rpos rneg], 'regioncenters', 'colormap');
    
    drawnow, snapnow
    
    r_medium = [rpos rneg];             % save with region names, also for data extraction below
    
    r_medium = extract_data_and_plot(r_medium, MC_Setup); % extract data from regions
    
else
    
    disp('No significant results.');
    
end

%% Lenient .05 cluster-extent corrected: Table and slice montage with autolabeled regions
% Excluding voxels already reported at more stringent/localizable thresholds

printhdr('Voxel-wise FWER + lenient primary threshold (.05) cluster-extent');

if ~isempty(r_lenient)
    
    %warning off, r = cluster2region(cl{4}); warning on
    
    montage(r_lenient, 'colormap');
    create_figure('surface'); axis off; surface(r_lenient, 'cutaway');
    
    drawnow, snapnow
    
    [rpos, rneg] = table(r_lenient);
    montage([rpos rneg], 'regioncenters', 'colormap');
    
    drawnow, snapnow
    
    r_lenient = [rpos rneg];             % save with region names, also for data extraction below
    
    r_lenient = extract_data_and_plot(r_lenient, MC_Setup); % extract data from regions
    
    
else
    
    disp('No significant results.');
    
end

close all       % save memory, etc.


% Save region_objects.mat file with regions and extracted data

save region_objects r_height r_stringent r_medium r_lenient


%% All thresholds: Montage and individual points
%
% % Reload masked activation map and display slice montage
% % The code below uses the CANlab object-oriented tools, in
% % CANlab_Core_Tools repository. There are many more options for
% % visualization and data analysis

disp('Activation_FWE_all.img')

% load saved results mask (union of all thresholds)
img = fmri_data('Activation_FWE_all.img', 'noverbose');  % Create an fmri_data-class object

create_figure('slice montage'); axis off
o2 = montage(img);                                       % Display montage, and return an fmridisplay object called o2

drawnow, snapnow

% Print transparent blobs and activation points on slice montage

% load database with coordinates
load SETUP DB
if ~isfield(DB, 'xyz'), DB.xyz = [DB.x DB.y DB.z]; end

o2 = removeblobs(o2);
o2 = addblobs(o2, r_lenient, 'trans', 'transvalue', .6);
o2 = addpoints(o2, DB.xyz, 'Color', 'none', 'MarkerFaceColor', [1 1 0], 'Marker', 'o', 'MarkerSize', 4);

drawnow, snapnow

o2 = removeblobs(o2);

drawnow, snapnow

close all



%% Surface renderings: Points on surfaces

% No longer needed
% Create a cutaway surface with blobs
% surface(img, 'cutaway', 'ycut_mm', -30);
% drawnow, snapnow

% Plot points as spheres on a canonical surface:

plot_points_on_surface2(DB.xyz, {[1 1 0]});

drawnow, snapnow

close all


%% Neurochemistry: PET binding map similarity

% [pet, petlabels] = load_image_set('receptorbinding');
% 
% stats = image_similarity_plot(h, 'mapset', pet, 'cosine_similarity', 'bicolor', 'networknames', pet.metadata_table.target, 'plotstyle', 'polar', 'average');
% 
% title('PET tracer profile');




% SUBFUNCTIONS
% ---------------------------------------------------------------------------


function region_obj = extract_data_and_plot(region_obj, MC_Setup)

indicator_maps = MC_Setup.unweighted_study_data;

if ~isempty(region_obj)
    
    % Extract indicator for contrasts that activate within 10 mm (5 vox) of
    % significant voxels region object r
    studybyroi = Meta_cluster_tools('getdata', region_obj, indicator_maps, MC_Setup.volInfo);
    
    for i = 1:length(region_obj)
        region_obj(i).dat = studybyroi(:, i);
    end
    
    % studybyroi: contrasts x regions, values 1/0 for whether each contrast activates the region
    % studybyset: contrasts x 1, values 1/0 for whether each contrast activates any region in the set
    
    create_figure('Activation_proportions')
    
    [n, k] = size(studybyroi);
    prop_by_condition = sum(studybyroi) ./ n;                        % Proportion of contrasts activating each ROI
    se = ( (prop_by_condition .* (1-prop_by_condition) ) ./ n ).^.5; % Standard Error based on binomial distribution
    
    han = bar(prop_by_condition);
    ehan = errorbar(prop_by_condition, se);
    
    labels = format_strings_for_legend({region_obj.shorttitle});
    
    set(han, 'EdgeColor', [0 0 .5], 'FaceColor', [.3 .3 .6]);
    set(ehan, 'Color', [0 0 .5], 'LineStyle', 'none', 'LineWidth', 3);
    set(gca, 'XTick', 1:k, 'XLim', [.5 k+.5], 'XTickLabel', labels);
    set(gca, 'XTickLabelRotation', 70);
    ylabel('Proportion of studies activating');
    
end % not empty

end % function
