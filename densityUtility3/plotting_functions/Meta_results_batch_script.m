%% CANlab MKDA meta-analysis results report
% "What's in the box": A series of tables and plots for a standard
% coordinate-based meta-analysis created with Meta_Activation_FWE in the
% CANlab MKDA meta-analysis toolbox.
% 
% The script Meta_results_batch_script is meant to be run stand-alone or
% used in publish_meta_analysis_report to create an HTML report with
% results and plots.

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

%% Voxel-wise FWER: Table and slice montage with autolabeled regions

if ~isempty(cl{1})
    
    warning off, r = cluster2region(cl{1}); warning on
    [rpos, rneg] = table(r);
    montage([rpos rneg], 'regioncenters', 'colormap');
    
    drawnow, snapnow
    
else
    
    disp('No significant results.');
    
end

[r_fwer] = [rpos rneg];             % save for data extraction below

%% Stringent .001 cluster-extent corrected: Table and slice montage with autolabeled regions

if ~isempty(cl{2})
    
    r = cluster2region(cl{2});
    [rpos, rneg] = table(r);
    montage([rpos rneg], 'regioncenters', 'colormap');
    
    drawnow, snapnow
    
else
    
    disp('No significant results.');
    
end

%% Stringent .001 cluster-extent corrected: Table and slice montage with autolabeled regions
% Excluding voxels already reported at more stringent/localizable thresholds

if ~isempty(cl{2})
    
    warning off, r = cluster2region(cl{2}); warning on
    [rpos, rneg] = table(r);
    montage([rpos rneg], 'regioncenters', 'colormap');
    
    drawnow, snapnow
    
else
    
    disp('No significant results.');
    
end

%% Medium .01 cluster-extent corrected: Table and slice montage with autolabeled regions
% Excluding voxels already reported at more stringent/localizable thresholds

if ~isempty(cl{3})
    
    warning off, r = cluster2region(cl{3}); warning on
    [rpos, rneg] = table(r);
    montage([rpos rneg], 'regioncenters', 'colormap');
    
    drawnow, snapnow
    
else
    
    disp('No significant results.');
    
end

%% Lenient .05 cluster-extent corrected: Table and slice montage with autolabeled regions
% Excluding voxels already reported at more stringent/localizable thresholds

if ~isempty(cl{4})
    
    warning off, r = cluster2region(cl{4}); warning on
    [rpos, rneg] = table(r);
    montage([rpos rneg], 'regioncenters', 'colormap');
    
    drawnow, snapnow
    
else
    
    disp('No significant results.');
    
end

close all       % save memory, etc.

%% Reload masked activation map and display slice montage
% The code below uses the CANlab object-oriented tools, in
% CANlab_Core_Tools repository. There are many more options for
% visualization and data analysis

% load saved results mask (union of all thresholds)
img = fmri_data('Activation_FWE_all.img', 'noverbose');  % Create an fmri_data-class object

create_figure('slice montage'); axis off
o2 = montage(img);                                       % Display montage, and return an fmridisplay object called o2

drawnow, snapnow

%% Print transparent blobs and activation points on slice montage

% load database with coordinates
load SETUP DB
if ~isfield(DB, 'xyz'), DB.xyz = [DB.x DB.y DB.z]; end

o2 = removeblobs(o2);   
o2 = addblobs(o2, r, 'trans', 'transvalue', .4);
o2 = addpoints(o2, DB.xyz, 'MarkerFaceColor', [.5 0 0], 'Marker', 'o', 'MarkerSize', 4);

drawnow, snapnow

close all

%% Surface renderings: Cutaways and points on surfaces 

% Create a cutaway surface with blobs
surface(img, 'cutaway', 'ycut_mm', -30);

drawnow, snapnow

% Plot points as spheres on a canonical surface:

plot_points_on_surface2(DB.xyz, {[.5 0 0]});

drawnow, snapnow

close all

%% Study proportion data for voxel-wise FWER-corrected regions
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
                       
indicator_maps = MC_Setup.unweighted_study_data;

if ~isempty(r_fwer)
    
    % Extract indicator for contrasts that activate within 10 mm (5 vox) of
    % significant voxels region object r
    [studybyroi,studybyset] = Meta_cluster_tools('getdata', r_fwer, indicator_maps, MC_Setup.volInfo);
    
    % studybyroi: contrasts x regions, values 1/0 for whether each contrast activates the region
    % studybyset: contrasts x 1, values 1/0 for whether each contrast activates any region in the set
    
    create_figure('Activation_proportions')
    
    [n, k] = size(studybyroi);
    prop_by_condition = sum(studybyroi) ./ n;                        % Proportion of contrasts activating each ROI
    se = ( (prop_by_condition .* (1-prop_by_condition) ) ./ n ).^.5; % Standard Error based on binomial distribution
    
    han = bar(prop_by_condition);
    ehan = errorbar(prop_by_condition, se);
    
    set(han, 'EdgeColor', [0 0 .5], 'FaceColor', [.3 .3 .6]);
    set(ehan, 'Color', [0 0 .5], 'LineStyle', 'none', 'LineWidth', 3);
    set(gca, 'XTick', 1:k, 'XLim', [.5 k+.5], 'XTickLabel', {r_fwer.shorttitle});
    set(gca, 'XTickLabelRotation', 70);
    ylabel('Proportion of studies activating');
    
end


