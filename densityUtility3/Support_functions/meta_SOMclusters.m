function varargout = meta_SOMclusters(meth,SOMResults,varargin)
    % varargout = meta_SOMclusters(meth,SOMResults,varargin)
    %
    % Multi-function toolbox for working with sets of clusters,
    % particularly those extracted from SOM parcellation
    %
    % [cl,anyStudy,studyByCluster] = meta_SOMclusters(SOMResults,theData,whSOM,meth,[Xi])
    %
    % Make clusters and indicator matrices for SOM results
    %
    % Methods:
    % -----------------------------------------------------------------------
    % 'clusters'
    % SOMResults = meta_SOMclusters('clusters',SOMResults,theData)
    %   Displays SOM Results on brain orthviews
    %   adds names, colors, and clusters (cl) fields to SOMResults
    %   Uses existing names and colors if already present
    %   Returns anyStudy field -> studies x SOMs, does any study activate in
    %   the SOM cluster?
    % -----------------------------------------------------------------------
    % 'single'
    % [cl,anyStudy,studyByCluster] = meta_SOMclusters('single',SOMResults,theData,whSOM)
    %   Gets anyStudy and studyByCluster for a single SOM cluster
    %
    % -----------------------------------------------------------------------
    % 'mca'
    % SOMResults = meta_SOMclusters('mca',SOMResults,[whdata])
    % Uses nometric MDS to get coordinate space
    % Uses permutation test testclust to get classes (superclusters)
    % Returns new region colors in SOMResults.MCA.colors
    % whdata can be 'som' (default) to run on SOM averages
    % or 'cluster' to run on individual cluster data (SOMs may have
    % multiple contiguous clusters).
    % -----------------------------------------------------------------------
    % 'stats'
    % SOMResults = meta_SOMclusters('stats',SOMResults,niterations,[whdata])
    % Uses permutation test to test for significant relationships among
    % regions.  Returns uncorrected and fdr-corrected p-values and
    % signficance matrices.
    % % whdata can be 'som' (default) to run on SOM averages
    % or 'cluster' to run on individual cluster data (SOMs may have
    % multiple contiguous clusters).
    % output returned in SOMResults.stats
    % -----------------------------------------------------------------------
    % 'names'
    % SOMResults = meta_SOMclusters('names',SOMResults,[addflag])
    %   Assign text names to SOMs
    %
    % -----------------------------------------------------------------------
    % 'xfig'
    % f1 = meta_SOMclusters('xfig',SOMResults,[colors])
    %   Make a graph (nmdsfig.m) of regions; no lines
    %   Uses SOMResults.colors if not entered separately
    % -----------------------------------------------------------------------
    % 'surface'
    % surfhandle = meta_SOMclusters('surface',SOMResults,keywd,radius)
    %   Make a surface plot of regions on any surface
    %   keywd can be any of the keywords used by addbrain
    %   i.e., empty for surface, 'left', 'right', 'thalamus','brainstem',
    %   etc.
    %   Uses SOMResults.MCA.colors if present, or SOMResults.colors
    %   otherwise
    %
    % ----------------------------------------------------------------------
    % 'somstats'
    % SOMResults = meta_SOMclusters('somstats',SOMResults,[whfield])
    % Needs fields: Xi, w, anyStudy (default; or other field of your choice)
    % Creates fields: prop_by_som, matrix of proportion of activations for SOM x task
    %               upper_ci, binomial upper confidence half-length
    %
    % To run for SOM-task associations, use:
    % SOMResults = meta_SOMclusters('somstats',SOMResults);
    %
    % To run for cluster-task associations, use:
    % SOMResults = meta_SOMclusters('somstats',SOMResults,'studyByCluster');
    %
    %
    % ----------------------------------------------------------------------
    % 'sombar'
    % [SOMResults,han] = meta_SOMclusters('sombar',SOMResults,[tasks,soms])
    % Needs fields: prop_by_som, upper_ci, alltasknms, names
    % Creates barplot of activation proportions for selected SOMs and
    % selected tasks
    % ----------------------------------------------------------------------
    % 'somcontrast'
    % meta_SOMclusters('somcontrast',SOMResults,[tasks])
    % Uses chi-square analysis to find SOMs that differ across conditions.
    % Show bar plot and brain slice plots of the significant SOMs.
    %
    % If no tasks are entered, gives you a menu to choose them.
    %
    % Data used:
    % Matrix sizes:     uses SOMResults.prop_by_som
    % Activation data:  default: .anyStudy; backup: .studyByCluster
    % Tasks:            .Xi
    % Weights:          .w
    %
    % To switch from using anyStudy to StudyByCluster, run somstats
    % subfunction (see above) with that field first.
    %
    % tor wager, Aug. 06

    cl = []; anyStudy = []; studyByCluster = [];

    if ~(exist('meth') == 1) || isempty(meth), meth = 'group';, end

    switch meth
        case 'clusters'

            theData = varargin{1};
            SOMResults = group_orthviews(SOMResults,theData);

            varargout{1} = SOMResults;


        case 'single'

            % [cl,anyStudy,studyByCluster] =
            % meta_SOMclusters(SOMResults,theData,whSOM,meth,Xi

            theData = varargin{1};
            whSOM = varargin{2};
            [cl,anyStudy,studyByCluster] = get_single_som(SOMResults,theData,whSOM);

            varargout{1} = cl;
            varargout{2} = anyStudy;
            varargout{3} = studyByCluster;


        case 'mca'

            if length(varargin) == 0, whdata = 'som'; else whdata = varargin{1}; end

            switch lower(whdata)
                case 'som'
                    whfield = 'anyStudy';
                case 'cluster'
                    whfield = 'studyByCluster';
                otherwise error('whdata must be som or cluster');
            end

            if ~isfield(SOMResults,whfield) || isempty(SOMResults.(whfield))
                error(['You must have ' whfield ' field in SOMResults.  Try ''clusters'' first.'])
            end

            fprintf(1,'Performing multiple correspondence (MCA) on clusters.\n');
            data = double(SOMResults.(whfield));
            SOMResults.MCA = tor_mca(data,'phi',[],[],'nmds'); % changed from 'correspondence', May 5, 2007; n flag changed from 1 to []
            varargout{1} = SOMResults;

        case 'stats'

            niter = varargin{1};
            whdata = 'som';
            if length(varargin) > 1, whdata = varargin{2}; end

            switch lower(whdata)
                case 'som'
                    whfield = 'anyStudy';
                case 'cluster'
                    whfield = 'studyByCluster';
                otherwise error('whdata must be som or cluster');
            end

            disp(['Getting data from ' whfield ])
            if ~isfield(SOMResults,whfield) || isempty(SOMResults.(whfield))
                error(['You must have ' whfield ' field in SOMResults.  Try ''clusters'' first.'])
            end

            if ~isfield(SOMResults,'Xi') || isempty(SOMResults.Xi)
                error('You must include task indic matrix Xi as a field in SOMResults.  Try ''Meta_Task_Indicator'' first.')
            end

            dat = full(double(SOMResults.(whfield)));
            nbrain = size(dat,2);
            dat = [dat SOMResults.Xi];
            SOMResults.stats = permute_mtx(dat,'ss',1,niter,nbrain);
            varargout{1} = SOMResults;



        case 'names'
            if length(varargin) > 0,addflag = varargin{1}; else, addflag = 0; end
            SOMResults.names = som_names(SOMResults.cl,addflag);
            varargout{1} = SOMResults;



        case 'xfig'
            %nmdsfig(MCA.score(:,1:2),ones(size(MCA.score,1),1),OUT.names,OUT.sigbonf);
            n = size(SOMResults.MCA.coords,1);

            if length(varargin) > 0
                disp('Using custom input colors.')
                colors = varargin{1};
            elseif isfield(SOMResults,'MCA') && isfield(SOMResults.MCA,'colors')
                disp('Using MCA colors.')
                colors = SOMResults.MCA.colors;
            else
                disp('Using initial random region colors.')
                colors = SOMResults.colors;
            end

            names = SOMResults.names;
            sig = zeros(n);
            if isfield(SOMResults,'stats') && size(SOMResults.stats.sigfdr,1) == n
                sig = SOMResults.stats.sigfdr;
            else
                disp(['SOMResults.stats is either missing or does not seem to match this data.'])
                disp(['Point sizes will not be based on stats'])
            end

            if length(names) == n
                f1 = nmdsfig(SOMResults.MCA.coords(:,1:2),'classes',(1:n)','names',SOMResults.names,'sig',sig, ...
                    'sizescale',[4 16],'nolines','colors',colors); axis off;
            else
                disp(['SOMResults.names is wrong length for this plot.  Omitting names.']);
                f1 = nmdsfig(SOMResults.MCA.coords(:,1:2),'classes',(1:n)','sig',sig, ...
                    'sizescale',[4 16],'nolines','colors',colors); axis off;
            end

            % this would use default colors
            %f1 = nmdsfig(SOMResults.MCA.coords(:,1:2),'classes',MCA.ClusterSolution.classes,'names',SOMResults.names,'sig',SOMResults.stats.sigfdr, ...
            %'sizescale',[4 16],'nolines'); axis off;

            varargout{1} = f1;


        case 'somstats'

            whfield = 'anyStudy';
            if length(varargin) > 0, whfield = varargin{1}; end
            SOMResults = som_task_stats(SOMResults,whfield);
            varargout{1} = SOMResults;

        case 'sombar'

            tasks = []; soms = [];
            if length(varargin) > 0, tasks = varargin{1}; end
            if length(varargin) > 1, soms = varargin{2}; end

            [SOMResults,han] = som_barplot(SOMResults,tasks,soms);
            varargout{1} = SOMResults;
            varargout{2} = han;

        case 'somcontrast'
            tasks = []; soms = [];
            if length(varargin) > 0, tasks = varargin{1}; end
            SOM_contrast(SOMResults,tasks);

        case 'surface'
            keywd = varargin{1};
            if length(varargin) > 1, r = varargin{2}; else r = 2; end
            sh = make_surface(SOMResults,keywd,r);
            varargout{1} = sh;

        otherwise
            error('Unknown method.  see help file for choices.');
    end


    return



% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% Subfunction

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------




function [cl,anyStudy,studyByCluster] = get_single_som(SOMResults,theData,whSOM)
    %whSOM = 1;

    hdr = SOMResults.header;
    vol = zeros(hdr.dim(1:3));
    ii = find(SOMResults.IDX == whSOM);
    vol(SOMResults.iMask(ii)) = SOMResults.WTS(ii);
    cl = mask2clusters(vol,hdr.mat);

    % view clusters
    %cluster_orthviews(cl,'unique');

    % figure out whether any study produced a non-zero value in the SOM area
    x = theData(ii,:);
    anyStudy = (sum(x,1) > 0)';

    if nargout > 2
        % figure out whether any study produced a non-zero value in each contiguous
        % cluster within the SOM
        [x,y,z]=ind2sub(hdr.dim(1:3),SOMResults.iMask(ii));
        indx = spm_clusters([x y z]');

        for i = 1:max(indx)
            wh = ii(indx == i); % which points in cluster
            x = theData(wh,:);  % pts by studies
            studyByCluster(:,i) = (sum(x,1) > 0)';
        end
        studyByCluster = double(studyByCluster);

        %burt = double(studyByCluster)' * double(studyByCluster);
    end

    return




% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% Subfunction

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


function SOMResults = group_orthviews(SOMResults,theData)
    % Displays SOM Results on brain orthviews
    % adds names, colors, and clusters (cl) fields to SOMResults

    fprintf(1,'meta_SOMclusters: Getting and plotting all SOMs. 000');
    cl = {};

    % names and sizes
    n = size(SOMResults.SOM,2);

    % names
    if ~isfield(SOMResults,'names') || isempty(SOMResults.names)
        for i=1:n
            SOMResults.names{i}=[num2str(i)];
        end
    end
    % colors
    if ~isfield(SOMResults,'colors') || isempty(SOMResults.colors)
        for i=1:n
            SOMResults.colors{i} = rand(1,3);
        end
    end

    for i = 1:n
        fprintf(1,'\b\b\b%03d',i);
        [cltmp,anyStudy(:,i),studyByCluster{i}] = get_single_som(SOMResults,theData,i);
        cltmp = rmfield(cltmp,'P');

        % remove small clusters (less than 20)
        wh = cat(1,cltmp.numVox) < 20;
        cltmp(wh) = [];
        studyByCluster{i}(:,wh) = [];

        % save soms by cluster
        sombycl{i} = repmat(i,1,length(cltmp));

        cl{i} = cltmp;
        nvox(i) = sum(cat(1,cl{i}.numVox));

        if i == 1, cluster_orthviews(cl{i},SOMResults.colors(i),'solid');
        else cluster_orthviews(cl{i},SOMResults.colors(i),'add','solid');
        end
    end

    SOMResults.cl = cl;
    SOMResults.anyStudy = anyStudy;

    SOMResults.studyByCluster = cat(2,studyByCluster{:});
    SOMResults.sombycl = cat(2,sombycl{:});

    return



% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% Subfunction

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


function names = som_names(cl,addflag)
    % function names = cluster_names(cl,[addflag])
    %
    % Assign names to cl(x).shorttitle
    % Do not use spaces or underscores or special chars for best results
    %
    % Addflag is optional; if 1, uses current orthview display

    spm_orthviews('Xhairs','on');


    for i = 1:length(cl)
        if addflag
            % don't create new figure
        else
            cluster_orthviews(cl{i},{[1 0 0]});
        end

        % show centers
        for j = 1:length(cl{i})
            spm_orthviews('Reposition',cl{i}(j).mm_center);
        end

        names{i} = input('Enter short name for this SOM set: ','s');

    end


    return



% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% Subfunction

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


function SOMResults = som_task_stats(SOMResults,whfield)
    % SOMResults = meta_SOMclusters('somstats',SOMResults)
    % Needs fields: Xi, w, anyStudy
    % Creates fields: prop_by_som, matrix of proportion of activations for SOM x task
    %               upper_ci, binomial upper confidence half-length

    if nargin < 2, whfield = 'anyStudy'; end

    nvars = size(SOMResults.(whfield),2);
    ntasks = size(SOMResults.Xi,2);

    % get stats for the entire matrix of SOMs
    [icon,ctxtxi,betas,num_by_condition,prop_by_condition] = meta_apply_contrast(full(double(SOMResults.(whfield)')), ...
        SOMResults.Xi,SOMResults.w,ones(1,ntasks));

    % hyp test: just an example, not valid
    %pvals = binomcdf(num_by_condition,repmat(sum(Xi),nvars,1),.5);

    % get confidence intervals based on estimated proportions
    % shrink towards .5 if proportion-hat is very extreme
    % ub is upper bound

    mytotals = repmat(sum(SOMResults.Xi),nvars,1);
    ub = binoinv(.975,mytotals,min(.9,max(prop_by_condition,.1)));
    ub = ub ./ mytotals;
    ub = ub - prop_by_condition;

    SOMResults.prop_by_som = prop_by_condition;
    SOMResults.upper_ci = ub;

    return



% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% Subfunction

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


function [SOMResults,han] = som_barplot(SOMResults,tasks,soms)

    if ~isfield(SOMResults,'prop_by_som')
        SOMResults = meta_SOMclusters('somstats',SOMResults);
    end

    if nargin < 2 || isempty(tasks)
        for i = 1:length(SOMResults.alltasknms)
            fprintf(1,'%3.0f\t%s\n',i,SOMResults.alltasknms{i});
        end
        tasks = input('Enter indices of tasks to plot: ');
    end

    if nargin < 3 || isempty(soms)
        for i = 1:length(SOMResults.names)
            fprintf(1,'%3.0f\t%s\n',i,SOMResults.names{i});
        end
        soms = input('Enter indices of SOMs to plot: ');
    end

    tasknames = SOMResults.alltasknms(tasks);

    nvars = size(SOMResults.prop_by_som,1);
    names = get_somnames(SOMResults,nvars);

    somnames = names(soms);

    myprops = SOMResults.prop_by_som(soms,tasks);
    myci = SOMResults.upper_ci(soms,tasks);

    tor_fig;
    han = barplot_grouped(myprops,myci,[],tasknames,'inputmeans');
    set(gca,'XTickLabel',somnames);
    ylabel('Proportion of studies activating');

    scn_export_papersetup(300)
    return

function names = get_somnames(SOMResults,nvars)
    % get names of vars we're interested in
    %if we're using clusters instead of SOMs,
    % length will not match; create names with numbers
    names = SOMResults.names;

    if length(names) ~= nvars
        disp(['SOMResults.names does not match.  Assuming you''re using clusters.'])
        for i = 1:nvars
            names{i} = num2str(i);
        end
    end
    return

    
    
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% Subfunction

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


function SOM_contrast(SOMResults,tasks)

    if nargin < 2 || isempty(tasks)
        for i = 1:length(SOMResults.alltasknms)
            fprintf(1,'%3.0f\t%s\n',i,SOMResults.alltasknms{i});
        end
        tasks = input('Enter indices of tasks to compare with chi-square analysis: ');
    end
    conname = input('Name this contrast (no spaces/special chars): ','s');

    
    % Matrix sizes:     uses SOMResults.prop_by_som
    % Activation data:  default: .anyStudy; backup: .studyByCluster
    % Tasks:            .Xi
    % Weights:          .w
    
    % names and sizes
    % ---------------------------------------------------------------------

    nvars = size(SOMResults.prop_by_som,1);
    ntasks = length(tasks);
    tasknames = SOMResults.alltasknms(tasks);
    somnames = get_somnames(SOMResults,nvars);

    % figure out which data field to use based on sizes of prop_by_som
    % get data and cell array of all clusters (based on cl or soms) to use
    % for display
    % ---------------------------------------------------------------------

    if nvars == size(SOMResults.SOM,2)
        whfield = 'anyStudy';
        disp(['Testing based on whether studies activated in SOMs.']);
        cl = SOMResults.cl;
    elseif nvars == size(SOMResults.studyByCluster,2)
        whfield = 'studyByCluster';
        disp(['Testing based on whether studies activated in individual clusters (parcels).']);
        % put in cell array form for compatibility
        cltmp = cat(2,SOMResults.cl{:});
        for i = 1:length(cltmp)
            cl{i} = cltmp(i);
        end
    else
        error('I can''t figure out where the data in prop_by_som came from.  It doesn''t match either clusters or SOMs.')
    end
    disp(['Getting data from ' whfield]);
    dat = double(full(SOMResults.(whfield)));


    % do Chi-sq analyses
    % ---------------------------------------------------------------------

    condf = SOMResults.Xi(:,tasks) .* repmat(1:ntasks,size(SOMResults.Xi,1),1);
    wh = sum(condf,2) - max(condf,[],2);
    if any(wh), warning('Task categories compared are not mutually exclusive.'), condf(wh,:) = 0; end
    condf = sum(condf,2);
    % eliminate empty rows
    wh = condf == 0;

    w = SOMResults.w;
    w = w./mean(w);
    w(wh) = [];
    sig = [];

    for som = 1:nvars
        y = [dat(:,som) condf];
        y(wh,:) = [];
        [chi2(som),df(som),chi2p(som),sig(som),warn(som),tab] = chi2test(y,'obs',w,1);

        % these should match what's in SOMResults.props_by_som, but we
        % don't need them because they do.
        %props(som,:) = tab(2,:) ./ sum(tab);

        if ntasks == 2 && size(tab,1) > 1
            % positive or negative; + means 1st col is greater, neg means
            % 2nd col is greater
            sig(som) = sig(som) .* -sign(diff(tab(2,:) ./ sum(tab)));
        end

    end

    % barplot
    % ---------------------------------------------------------------------
    sigsoms = find(sig);
    if isempty(sigsoms), disp('No significant results.'), return, end
    
    [SOMResults,han] = meta_SOMclusters('sombar',SOMResults,tasks,sigsoms);
    name = [conname '_barplot'];
    scn_export_papersetup(300);
    saveas(gcf,name,'png');

    % show clusters
    % ---------------------------------------------------------------------
    sigsoms = find(sig>0);
    addstr = 'noadd';
    if ~isempty(sigsoms)
        mycl = cat(2,cl{sigsoms});
        cluster_orthviews(mycl,{[0 0 0]},'solid');
        addstr = 'add';
    end

    sigsoms = find(sig<0);
    if ~isempty(sigsoms)
        mycl = cat(2,cl{sigsoms});
        cluster_orthviews(mycl,{[1 1 1]},addstr,'solid');
    end

    sigsoms = find(sig);
    mycl = cat(2,cl{sigsoms});
    %sigcl = mycl(sigsoms);

    
    % tables
    dotables = input('Print tables? (1 or 0) ');
    if dotables
        prop = SOMResults.prop_by_som(:,tasks) ;

        % header
        fprintf(1,'SOMcl\tX\tY\tZ\tVox\t');
        for j = 1:ntasks, fprintf(1,'%s\t',tasknames{j}); end
        fprintf(1,'Chi2\tp\t\n');
            
        % body
        mycl = cat(2,SOMResults.cl{:});
        ncl = length(cl);
        for i = 1:ncl
            if chi2p(i) < .005, sigstr = '***';
            elseif chi2p(i) < .01; sigstr = '**';
            elseif chi2p(i) < .05, sigstr = '*';
            else, sigstr = '';
            end
            
            if ~isfield(mycl(i),'shorttitle') || isempty(mycl(i).shorttitle), mycl(i).shorttitle = ['Cl ' num2str(i)]; end
            
            fprintf(1,'%s\t%3.0f\t%3.0f\t%3.0f\t%3.0f\t',mycl(i).shorttitle,mycl(i).mm_center(1),mycl(i).mm_center(2),mycl(i).mm_center(3),mycl(i).numVox);
             for j = 1:ntasks, fprintf(1,'%3.0f\t',100.*prop(i,j)); end
             fprintf(1,'%3.2f\t%3.4f%s\t\n',chi2(i),chi2p(i),sigstr);
        end
        
    end


    

    dofigs = input('Save brain slices and surfaces? (1 or 0) ');
    if ~isempty(sigsoms) && dofigs

        % get list of color indices (1 or 2) for each cluster
        % we need to do this for the SOM option only
        if nvars == size(SOMResults.SOM,2)
            sig2 = [];
            for i = 1:length(sigsoms)
                len = length(cl{sigsoms(i)});
                sig2 = [sig2 repmat(sig(sigsoms(i)),1,len)];
            end
            sig = sig2;
        end

        % surface: lateral
        colorindx = (sig(find(sig)) > 0) + 1;
        colors = {[1 1 1] [0 0 0]};     % neg then pos
        tor_fig; sh = cluster_surf(mycl(1),colors(colorindx(1)));
        for i = 2:length(mycl(sigsoms))
            cluster_surf(mycl(sigsoms(i)),2,colors(colorindx(i)),sh);
        end

        name = [conname '_lat_surf'];
        saveas(gcf,name,'fig');
        saveas(gcf,name,'png');

        % surface: brainstem
        tor_fig; sh = cluster_surf(mycl(1),2,colors(colorindx(1)),'brainstem');
        for i = 2:length(mycl)
            cluster_surf(mycl(i),colors(colorindx(i)),sh);
        end

        name = [conname '_brainstem_surf'];
        saveas(gcf,name,'fig');
        saveas(gcf,name,'png');

        % surface: left
        tor_fig; sh = cluster_surf(mycl(1),2,colors(colorindx(1)),'left');
        for i = 2:length(mycl)
            cluster_surf(mycl(i),colors(colorindx(i)),sh);
        end

        name = [conname '_left_surf'];
        saveas(gcf,name,'fig');
        saveas(gcf,name,'png');

        name = [conname '_left_lat'];
        view(270,5);
        lighting gouraud; lightRestoreSingle(gca)
        saveas(gcf,name,'png');

        % surface: right
        tor_fig; sh = cluster_surf(mycl(1),2,colors(colorindx(1)),'right');
        for i = 2:length(mycl)
            cluster_surf(mycl(i),colors(colorindx(i)),sh);
        end

        name = [conname '_right_surf'];
        saveas(gcf,name,'fig');
        saveas(gcf,name,'png');


        name = [conname '_right_lat'];
        view(90,5);
        lighting gouraud; lightRestoreSingle(gca)
        saveas(gcf,name,'png');


        overlay = which('scalped_single_subj_T1.img');

        % slices
        name = [conname '_axial'];
        cluster_orthviews_showcenters(mycl,'axial',overlay,0);
        scn_export_papersetup(800);
        saveas(gcf,name,'png')

        name = [conname '_saggital'];
        cluster_orthviews_showcenters(mycl,'saggital',overlay,0);
        scn_export_papersetup(800);
        saveas(gcf,name,'png')

        name = [conname '_coronal'];
        cluster_orthviews_showcenters(mycl,'coronal',overlay,0);
        scn_export_papersetup(800);
        saveas(gcf,name,'png')
    end


    %     % get largest cluster in each significant SOM and image it
    %     clplot = [];
    %     for i = 1:length(sigsoms)
    %         mycl = SOMResults.cl{sigsoms(i)};
    %         vs = cat(1,mycl.numVox);
    %         wh = vs == max(vs);
    %         if isempty(clplot), clplot = mycl(wh); else clplot = merge_clusters(clplot,mycl(wh)); end
    %     end

    return



% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% Subfunction

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


function sh = make_surface(SOMResults,keywd,radmm)

    if isempty(keywd)
        sh = addbrain;
    else
        sh = addbrain(keywd);
    end

    colors = SOMResults.colors;
    if isfield(SOMResults,'MCA')
        colors = SOMResults.MCA.colors;
        disp('Using MCA colors')
    else
        disp('Using SOM random colors');
    end

    set(sh,'FaceAlpha',1)
    lightRestoreSingle(gca)
    %sh = cluster_surf(SOMResults.cl{1},SOMResults.colors{1},2);
    for i = 1:length(SOMResults.cl)
        cluster_surf(SOMResults.cl{i},colors(i),radmm,sh);
    end

    axis image; lighting gouraud; axis vis3d
    lightRestoreSingle(gca)

    return

