function Meta_Multivariate(MC_Setup,DB,imgname)

% -----------------------------------------------------------------------

% Prepare data matrix
% Needs: MC_Setup.unweighted_study_data, 'Activation_FWE_all.img' in
% current directory

% -----------------------------------------------------------------------
disp('Preparing data.')
dat = MC_Setup.unweighted_study_data;

% Get only significant voxels
volInfo = iimg_read_img('Activation_FWE_all.img',2);
is_significant = volInfo.image_indx(MC_Setup.volInfo.wh_inmask);

dat = dat(is_significant,:);

% -----------------------------------------------------------------------

% Prepare task indicator matrix
% Needs: DB structure, data matrix

% -----------------------------------------------------------------------
disp('Getting task proportions.')
% Get task averages: Proportion of studies in each condition
% Voxels x task conditions
[Xi,alltasknms,condf,allti,testfield,prop_by_condition,num_by_condition] = Meta_Task_Indicator(DB,dat);

ncontig = max(volInfo.cluster);

%classes = clusterdata(prop_by_condition,2,'linkage','average');

nSOM = 25; nIter = 5000;

% Calculate the SOM

SOMResults        = SOM_CalculateMap(dat,nSOM,nIter);

% Store the header information - needed for writing out results as images.

SOMResults.header = volInfo;
SOMResults.iMask = volInfo.wh_inmask;

% Organize the data into super clusters.
%[SOMResults.SuperCluster SOMResults.nCluster] = SOM_SuperCluster(SOMResults.SOM);



% DO MCA (multidimensional scaling) on regions to get classes (superclusters)
SOMResults = meta_SOMclusters('mca',SOMResults);

% get colors baesd on classes (superclusters)
% done automatically in tor_mca
SOMResults.MCA.colors = nmdsfig_tools('class_colors',MCA.ClusterSolution.classes);

% make a network graph with new colors
f1 = meta_SOMclusters('xfig',SOMResults);

% show orthviews in new colors
SOMResults = meta_SOMclusters('clusters',SOMResults,dat);