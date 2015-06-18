function SOMResults = Meta_Multivariate(MC_Setup,DB,nSOM,nIter,varargin)
    % SOMResults = Meta_Multivariate(MC_Setup,DB,nSOM,nIter,[additional args])
    %
    % Additional args:
    % 'mask', then mask image name
    % 'data', followed by data matrix
    %
    % This is a clearing-house function to run multivariate-based analysis of
    % meta-analysis data (peak coordinates).  The main steps are:
    %
    % 1) parcellation of the brain into a number of contiguous regions (clusters) using a
    % self-organizing map algorithm (by Scott Peltier and Robert Welsh).  An
    % "SOM" is defined here as the set of clusters that belong to the same
    % canonical function (the same map).
    %
    % 2) MCA: multiple correspondence analysis, or a variant thereof, on the
    % SOMs (or, alternatively, on clusters).  This involves nonmetric MDS on
    % the associations between SOMs, clustering of the SOMs, and permutation
    % testing to determine the best number of "superclusters" and whether there
    % is significant clustering of SOMs.
    %
    % 3) stats: saving statistics on significant relationships between SOMs ([default]; or
    % clusters) based on permutation test.
    %
    % 4) some tables and output
    %
    % To prepare the input, you need to run:
    % a) DB = read_database
    % b) DB = Meta_Setup(DB); to get DB structure with point locations and
    % fields
    % c) Meta_Activation_FWE('Setup',DB)  to prepare the MC_Setup structure with
    % the data (with convolution depending on DB.radius)
    % d) Meta_Activation_FWE('MC....  and 'Results' to get
    % Activation_FWE_all.img file with significant regions to parcellate
    % (Or change to add your own mask)
    %
    % A note on processing the output:
    % -------------------------------------------------------------------------
    % meta_SOMclusters is a utility function that works with the results of the
    % SOMResults structure, and has many useful methods.
    % e.g.,
    % [SOMResults,han] = meta_SOMclusters('sombar',SOMResults,tasks,soms);
    % Makes a bar graph of each task of interest for a set of given SOMs of interest
    %
    % meta_SOMclusters('contrast',SOMResults)
    % lets you specify a contrast across task conditions and returns images of
    % SOMs or clusters showing significant differences among the tasks.
    %
    % Many other useful functions exist also, such as plot_points_on_brain,
    % plot_points_on_subctx, and others

    maskimg = 'Activation_FWE_all.img';
    for i = 1:length(varargin)
        if isstr(varargin{i})
            switch varargin{i}

                % functional commands
                case 'mask', maskimg = varargin{i+1};
                case 'data', dat = varargin{i+1};

                otherwise, warning(['Unknown input string option:' varargin{i}]);
            end
        end
    end

    % -----------------------------------------------------------------------

    % Prepare data matrix
    % Needs: MC_Setup.unweighted_study_data, 'Activation_FWE_all.img' in
    % current directory

    % -----------------------------------------------------------------------

    disp('Preparing data.')
    if ~exist('dat','var'), dat = MC_Setup.unweighted_study_data; end

    % Get only significant voxels
    volInfo = iimg_read_img(maskimg,2);
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

    if isfield(DB,'SubjectiveWeights')
        w = DB.rootn .* DB.SubjectiveWeights(DB.pointind);
    else
        w = DB.rootn;
    end
    ncontig = max(volInfo.cluster);

    disp(['Saving Meta_Multiv_Task_Data_Setup.mat'])
    save Meta_Multiv_Task_Data_Setup


    % -----------------------------------------------------------------------

    % Parcel using self-organizing map (SOM)
    % Needs: data matrix, nSOM, nIter, volInfo (from mask, internal)

    % -----------------------------------------------------------------------

    % Calculate the SOM
    SOMResults        = SOM_CalculateMap(prop_by_condition,nSOM,nIter);

    % Store the header information - needed for writing out results as images.
    SOMResults.header = volInfo;
    SOMResults.iMask = volInfo.wh_inmask;

    % Organize the data into super clusters.
    %[SOMResults.SuperCluster SOMResults.nCluster] = SOM_SuperCluster(SOMResults.SOM);

    % add data
    SOMResults.dat = dat;

    % add task info
    SOMResults.Xi = Xi;
    SOMResults.alltasknms = alltasknms;
    SOMResults.prop_by_condition = prop_by_condition;
    SOMResults.num_by_condition = num_by_condition;
    SOMResults.w = w;

    save SOMResults SOMResults
    disp(['Saving SOMResults.mat'])

    % show orthviews and get SOM x task indicator matrix needed for next steps
    SOMResults = meta_SOMclusters('clusters',SOMResults,dat);



    % -----------------------------------------------------------------------

    % Do MCA and make graph
    % Needs: SOMResults

    % -----------------------------------------------------------------------
    % DO MCA (multidimensional scaling) on regions to get classes (superclusters)
    SOMResults = meta_SOMclusters('mca',SOMResults);

    scn_export_papersetup(400);
    saveas(gcf,'testclust','png');
    close;

    scn_export_papersetup(400);
    saveas(gcf,'shepardplot','png');
    close;


    % get colors baesd on classes (superclusters)
    % done automatically in tor_mca
    SOMResults.MCA.colors = nmdsfig_tools('class_colors',SOMResults.MCA.ClusterSolution.classes);

    disp(['Saving SOMResults.mat with MCA field'])
    save SOMResults SOMResults

    % -----------------------------------------------------------------------

    % Get stats for correlations among regions, tasks, and region-task
    % relationships, based on permutation
    % Needs: SOMResults with anyStudy and Xi (task) fields, and n. iterations

    % -----------------------------------------------------------------------

    SOMResults = meta_SOMclusters('stats',SOMResults,nIter);

    % Get univariate stats for the relationship between SOM and each task
    % (fdr p-values are given by the nonparametric analysis; this is mostly for
    % plotting proportions and confidence intervals).
    SOMResults = meta_SOMclusters('somstats',SOMResults);

    disp(['Saving SOMResults.mat with stats and prop_by_som fields'])
    save SOMResults SOMResults

    % -----------------------------------------------------------------------

    % Output graphs

    % -----------------------------------------------------------------------

    % make a network graph with new colors
    % need MCA and stats
    f1 = meta_SOMclusters('xfig',SOMResults);

    scn_export_papersetup(500);
    saveas(f1,'nmdsfig','png')
    saveas(f1,'nmdsfig','fig')

    % show orthviews in new colors
    SOMResults = meta_SOMclusters('clusters',SOMResults,dat);

    pruneval = input('Prune clusters based on size?  Enter return or voxel count threshold: '); 
    if pruneval
        % prune clusters based on size
        for i = 1:length(SOMResults.cl)
            wh = cat(1,SOMResults.cl{i}.numVox);
            wh = find(wh < pruneval);
            if length(wh) == length(SOMResults.cl{i}), wh = []; end
            SOMResults.cl{i}(wh) = [];
        end

    end


    % table of output: which Classes/regions are associated with which tasks
    diary class_by_task_table.txt
    classes = SOMResults.MCA.ClusterSolution.classes;
    nclasses = length(unique(classes));
    for i = 1:nclasses
        wh = find(classes == i);
        disp(['Class ' num2str(i)]);
        print_matrix(SOMResults.stats.part12.sigfdr(wh,:),SOMResults.alltasknms,SOMResults.names(wh));
        disp('')
    end
    diary off

    overlay = which('spm2_single_subj_T1_scalped.img');
    
    % orthviews of regions significantly related to each task condition
    ntasks = size(SOMResults.Xi,2);
    for i = 1:ntasks
        wh = find(SOMResults.stats.part12.sigfdr(:,i));
        for j = 1:length(wh)
            if j==1,addstr='noadd'; else addstr ='add'; end
            cluster_orthviews(SOMResults.cl{wh(j)},SOMResults.MCA.colors(wh(j)),addstr,'solid');
        end

        mycl = cat(2,SOMResults.cl{wh});

        name = [SOMResults.alltasknms{i} '_axial'];
        cluster_orthviews_showcenters(mycl,'axial',overlay);
        scn_export_papersetup(800);
        saveas(gcf,name,'png')

        name = [SOMResults.alltasknms{i} '_saggital'];
        cluster_orthviews_showcenters(mycl,'saggital',overlay);
        scn_export_papersetup(800);
        saveas(gcf,name,'png')

        name = [SOMResults.alltasknms{i} '_coronal'];
        cluster_orthviews_showcenters(mycl,'coronal',overlay);
        scn_export_papersetup(800);
        saveas(gcf,name,'png')
    end


    return



