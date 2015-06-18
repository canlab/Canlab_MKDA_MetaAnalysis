function varargout = Meta_cluster_tools(meth,varargin)
% varargout = Meta_cluster_tools(meth,varargin)
%
% This function contains multiple tools for working with clusters
% derived from Meta_Activation_FWE and Meta_SOM tools
%
%
% ------------------------------------------------------
% extract data and print a table for a set of meta-analysis clusters
%
% cl = make_table(cl,MC_Setup,['plot'],['successive'])
% Example:
% load Valence_Neg-Pos_Pos_clusters; load MC_Info
%
% cl1 = Meta_cluster_tools('make_table',cl{1},MC_Setup);
% cl = Meta_cluster_tools('make_table',cl,MC_Setup,'successive');
%
% ------------------------------------------------------
%
% Print tables of which contrasts in DB database activate within clusters
% with full contrast/database information:
% Meta_cluster_tools('activation_table', DB, MC_Setup, cl, [testfield1], [testfield2]);
% Meta_cluster_tools('activation_table', DB, MC_Setup, cl(1), 'Modality2', 'PosNeg');
%
% DB should be database, after Meta_Setup.m has been run.
% MC_Setup should be structure from Meta_Activation_FWE.m with setup info
% cl is clusters structure of results, e.g., from chi-square, etc.
% testfield1 and 2 are optional names of fields to calculate percentage
% activation for
% a contingency table is made if two fields are entered.
%
% ------------------------------------------------------
% get data for studies within rois
% studybyroi is studies activating in each cluster in cl. operator is "any" voxel in cluster counts
% [studybyroi,studybyset] = Meta_cluster_tools('getdata',cl,dat,[volInfo])
% [studybyroi,studybyset] = Meta_cluster_tools('getdata',cl,MC_Setup.unweighted_study_data,MC_Setup.volInfo)
% ------------------------------------------------------
%
% ------------------------------------------------------
% count studies by condition and plot [optional]
%
% [prop_by_condition,se,num_by_condition,n] = Meta_cluster_tools('count_by_condition',dat,Xi,w,doplot,[xnames],[seriesnames], [colors])
%
% [prop_by_condition,se,num_by_condition,n] = Meta_cluster_tools('count_by_condition',studybyset,MC_Setup.Xi,MC_Setup.wts,1)
%
% [prop_by_condition,se,num_by_condition,n] = ...
%   Meta_cluster_tools('count_by_condition',studybyroi,MC_Setup.Xi,MC_Setup.wts,1, ...
%   {'Right' 'Left'},MC_Setup.connames(1:5),{[1 0 0] [0 1 0] [1 0 1] [1 1 0] [0 0 1]});
%
% Xi = SOMResults.Xi(:,9:13);
% nms = SOMResults.alltasknms(9:13)
% w = SOMResults.w;
% colors = {[1 0 0] [0 1 0] [1 0 1] [1 1 0] [0 0 1]};
% [prop_by_condition,se,num_by_condition,n] = ...
%   Meta_cluster_tools('count_by_condition',studybyroi,Xi,w,1, ...
%   {'Right' 'Left'},nms,colors);
%
% ------------------------------------------------------
%
% Example:
% ------------------------------------------------------
% Run an analysis with Meta_Chisq_new, and then use these tools to get
% plots of regions. The lines below run the entire analysis.
% R = Meta_Chisq_new('compute',MC_Setup,'mask',mask);
% R = Meta_Chisq_new('write',R);
% [studybyroi,studybyset] = Meta_cluster_tools('getdata',cl,R.dat,R.volInfo);
% [prop_by_condition,se,num_by_condition,n] = Meta_cluster_tools('count_by_condition',studybyset,R.Xi,R.w,1);
% ------------------------------------------------------

switch meth
    
    case 'make_table'
        
        % cl = make_table(cl,MC_Setup,[doplot],[successiveflag])
        %
        doplot = 0; dosuccessive = 0;
        if any(strcmp(varargin,'successive')), dosuccessive = 1; end
        if any(strcmp(varargin,'plot')), doplot = 1; end
        
        varargout{1} = make_table(varargin{1},varargin{2},doplot,dosuccessive);
        
    case 'getdata'
        
        %[studybyroi,studybyset] = Meta_cluster_tools('getdata',cl,MC_Setup.unweighted_study_data)
        %[studybyroi,studybyset] = getdata(cl,inputdata)
        if length(varargin) < 3
            [varargout{1},varargout{2},varargout{3}] = getdata(varargin{1},varargin{2});
        else
            [varargout{1},varargout{2},varargout{3}] = getdata(varargin{1},varargin{2},varargin{3});
        end
        
    case 'count_by_condition'
        %[prop_by_condition,se,num_by_condition,n] = Meta_cluster_tools('count_by_condition',dat,Xi,w,doplot,[xnames],[seriesnames], [colors])
        %[prop_by_condition,se,num_by_condition,n] = count_by_condition(dat,Xi,w,doplot)
        
        if length(varargin) < 4, varargin{4} = 0; end
        if length(varargin) < 5, varargin{5} = []; end
        if length(varargin) < 6, varargin{6} = []; end
        if length(varargin) < 7, varargin{7} = []; end  %colors
        [varargout{1},varargout{2},varargout{3},varargout{4}] = count_by_condition(varargin{1},varargin{2},varargin{3},varargin{4},varargin{5},varargin{6},varargin{7});
        
        
    case 'activation_table'
        activation_table(varargin{:});
        
    otherwise
        disp('unknown method string.  doing nothing.');
        
end


end



% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% Get data within rois

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


function cl = make_table(cl,MC_Setup,doplot,dosuccessive)

% uses Xi and Xinms to get differences among conditions
% if there is no Xi field, we don't have differences among
% conditions, so create a dummy one to get overall proportion
if ~isfield(MC_Setup,'Xi')
    MC_Setup.Xi = ones(size(MC_Setup.wts));
    MC_Setup.Xinms = {'Act'};
end

if dosuccessive
    for i = 1:length(cl)
        disp(['Cluster cell ' num2str(i)])
        cl{i} = get_props_subfcn(cl{i},MC_Setup,doplot);
    end
else
    cl = get_props_subfcn(cl,MC_Setup,doplot);
end


% build table function call
if dosuccessive
    estr = 'cl = cluster_table_successive_threshold(cl,5';
else
    estr = 'cluster_table(cl,1,0';
end

fnames = MC_Setup.Xinms;
for i = 1:length(fnames), fnames{i} = [fnames{i} '_prop']; end
nconds = length(fnames);

for i = 1:nconds
    estr = [estr ',''' fnames{i} ''''];
end
estr = [estr ');'];

% run table
eval(estr)

end

% dependent on above:
function cl = get_props_subfcn(cl,MC_Setup,doplot)
disp(['getting clusters for local maxima at least 10 mm apart']);
cl = subclusters_from_local_max(cl,10);
cl = merge_nearby_clusters(cl,10,'recursive');

% get proportion of points activating in each condition in each region
disp('Getting studies that activated in each region.')
[studybyroi,studybyset] = Meta_cluster_tools('getdata',cl,MC_Setup.unweighted_study_data,MC_Setup.volInfo);

disp('Counting studies by condition')
[prop_by_condition,se,num_by_condition,n] = Meta_cluster_tools('count_by_condition',studybyroi,MC_Setup.Xi,MC_Setup.wts,doplot);

% get field names for conditions
fnames = MC_Setup.Xinms;
for i = 1:length(fnames), fnames{i} = [fnames{i} '_prop']; end
nconds = length(fnames);

fprintf(1,'Adding field to cl: %s\n',fnames{:});

% store proportions in clusters for table printout and posterity
for i = 1:length(cl)
    for j = 1:nconds
        cl(i).(fnames{j}) = 100 * prop_by_condition(i,j);
    end
end

end


function [studybyroi,studybyset, cl] = getdata(cl,inputdata,varargin)

if length(varargin) > 0
    volInfo = varargin{1};
    %maskname = volInfo.fname;
else
    disp('Using default mask. if your data has a different set of voxels, enter volInfo as input.')
    maskname = which('scalped_avg152T1_graymatter_smoothed.img');
    volInfo = iimg_read_img(maskname);
end

n_inmask_in = size(inputdata, 1);
if n_inmask_in ~= volInfo.n_inmask
    fprintf('*****************************\nWARNING\n*****************************\n')
    fprintf('Voxels in input data set: %3.0f\nVoxels in volInfo: %3.0f\n', n_inmask_in, volInfo.n_inmask);
    fprintf('These must match!\n*****************************\n');
end

nrois = length(cl);
nstudies = size(inputdata,2);

studybyroi = false(nstudies,nrois);

for i = 1:nrois
    [imgvec,maskvec] = iimg_clusters2indx(cl(i),volInfo);    %maskname);
    
    dat = inputdata(maskvec,:);
    
    cl(i).Z = sum(full(dat'));
    cl(i).Z_descrip = 'Unweighted sum of activating studies';
    
    studybyroi(:, i) = any(dat,1)';
    
    cl(i).activating_comparisons = studybyroi(:, i);
    
end

[imgvec,maskvec] = iimg_clusters2indx(cl,volInfo);   %maskname);
dat = inputdata(maskvec,:);
studybyset = full(any(dat)');

studybyroi = full(studybyroi);

end


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% Count studies in each region by condition and plot if asked for

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

function [prop_by_condition,se,num_by_condition,n] = count_by_condition(dat,Xi,w,doplot,varargin)

if nargin < 4, doplot = 0; end

[nstudies,ntasks] = size(Xi);
if size(dat,1) ~= nstudies, dat = dat'; end
if size(dat,1) ~= nstudies, error('data size does not match Xi'); end
nregions = size(dat,2);

% get stats for the entire matrix of SOMs
[icon,ctxtxi,betas,num_by_condition,prop_by_condition] = meta_apply_contrast(dat', ...
    Xi,w,ones(1,ntasks));

n = sum(Xi);
n = repmat(n,nregions,1);
se = ( (prop_by_condition .* (1-prop_by_condition) ) ./ n ).^.5;

if doplot
    xnames = [];
    seriesnames = [];
    mycolors = [];
    if length(varargin) > 0, xnames = varargin{1}; end
    if length(varargin) > 1, seriesnames = varargin{2}; end
    if length(varargin) > 2, mycolors = varargin{3}; end
    
    create_figure('barplot');
    fprintf('Sample size is %3.0f\n', size(Xi, 1));
    
    barplot_grouped(prop_by_condition, se, xnames, seriesnames, 'inputmeans', 'colors', mycolors);
    
    ylabel('Proportion of studies activating');
end

end








function activation_table(DB, MC_Setup, cl, testfield, testfield2)
% Table header

%cl = database2clusters(DB, cl, DB.radius_mm);

% Get list of activating studies
[studybyroi, studybyset] = Meta_cluster_tools('getdata', cl, MC_Setup.unweighted_study_data,MC_Setup.volInfo);

% Get which fields are valid to use
N = fieldnames(DB);

for i = 1:length(N)
    if length(size(DB.(N{i}))) == 2 && all(size(DB.(N{i})) == size(DB.x))
        include(i, 1) = true;
    else
        include(i, 1) = false;
    end
end

N = N(include);
myz = '=======================================================';

for r = 1:length(cl)
    
    if ~isfield(cl(r), 'shorttitle') || isempty(cl(r).shorttitle)
        cl(r).shorttitle = sprintf('%3.0f', r);
    end
    
    fprintf('%s\nRegion %s\n%s\n%s\n', myz, cl(r).shorttitle, cl(r).title, myz);
    
    % Header
    for f = 1:length(N)
        fprintf('%s\t', N{f});
    end
    fprintf('\n');
    
    whcons = DB.pointind(logical(studybyroi(:, 1)));
    n = length(whcons);
    
    for s = 1:n
        
        for f = 1:length(N)
            
            myval = DB.(N{f})(whcons(s));
            
            if iscell(myval)
                fprintf('%s\t', myval{1});
            elseif myval == round(myval)
                fprintf('%3.0f\t', myval);
            else
                fprintf('%3.2f\t', myval);
            end
            
            %             % save data for contingencies
            %             if iscell(myval)
            %             mydata{s, f} = myval{1};
            %             else
            %                  mydata{s, f} = myval;
            %             end
            
        end  % f = field
        
        fprintf('\n');
        
    end % s = study/contrast
    
    fprintf('\n_________________________________________________________\n');
    
    
    % now contingency table if we have it
    if exist('testfield', 'var')
        
        mydata = DB.(testfield)(DB.pointind);
        names = unique(mydata');
        
        actcons = logical(studybyroi(:, 1));
        
        fprintf('%s\t%s\t%s\t%s\n', 'Condition', 'Total Cons', 'ACtive Cons', '% Active');
        
        for c = 1:length(names)
            mytotal = sum(strcmp(mydata, names{c}));
            myactive = sum(strcmp(mydata(actcons), names{c}));
            
            fprintf('%s\t%3.0f\t%3.0f\t%3.0f%%\n', names{c}, mytotal, myactive, 100*myactive/mytotal);
            
        end
        
    end
    
    if exist('testfield2', 'var')
        
        fprintf('\n');
        
        mydata = DB.(testfield2)(DB.pointind);
        names = unique(mydata');
        
        actcons = logical(studybyroi(:, 1));
        
        fprintf('%s\t%s\t%s\t%s\n', 'Condition', 'Total Cons', 'Active Cons', '% Active');
        
        for c = 1:length(names)
            mytotal = sum(strcmp(mydata, names{c}));
            myactive = sum(strcmp(mydata(actcons), names{c}));
            
            fprintf('%s\t%3.0f\t%3.0f\t%3.0f%%\n', names{c}, mytotal, myactive, 100*myactive/mytotal);
            
        end
        
        
        % two-way table of percentages
        fprintf('\nContingency table of proportions of activating studies in combined categories\n');
        
        mydata = DB.(testfield)(DB.pointind);
        
        mydata2 = DB.(testfield2)(DB.pointind);
        
        [indx, names] = string2indicator(mydata);
        [indx2, names2] = string2indicator(mydata2);
        mytotal = indx' * indx2; % total contrasts in each combo
        
        clear indx
        for j = 1:length(names)
            indx(:, j) = double(strcmp(mydata(actcons), names{j}));
        end
        
        clear indx2
        for j = 1:length(names2)
            indx2(:, j) = double(strcmp(mydata2(actcons), names2{j}));
        end
        
        myactive = indx' * indx2; % total activating contrasts in each combo
        
        myperc = myactive ./ mytotal;
        
        print_matrix(myperc, names2, names)
        
        
        create_figure('contingency');
        imagesc(100 .* myperc)
        colorbar
        cm = colormap_tor([.3 0 .7], [1 1 0]);
        colormap(cm)
        title('Percentage of studies activating by category');
        set(gca, 'XTick', 1:length(names2), 'XTickLabel', names2);
        set(gca, 'YTick', 1:length(names), 'YTickLabel', names);
        snapnow
        
    end
    
    
    
end % r = cl

end % function





function [N,con,testfield] = check_fields(DB,testfield)
if isfield(DB,'N'), N = DB.N;, elseif isfield(DB,'Subjects'), N = DB.Subjects; else, N = NaN*zeros(size(DB.x));, end
if isfield(DB,'Contrast'), con = DB.Contrast;,
else, con = NaN*zeros(size(DB.x));,
    disp('Warning!  You must have a field called DB.Contrasts for the table function to work properly.');
end

if ~isfield(DB,'connumbers'),
    error('No DB.connumbers field, which is required.  See Meta_Setup to create this field and set up analysis.');
end

% Define testfield (field to display in table)
if isempty(testfield), try load SETUP testfield, catch, testfield = input('Cannot load testfield from SETUP.  Type name of field in DB to display: ','s'), go = 1;, end, end

if isempty(testfield), testfield = input('Cannot load testfield from SETUP.  Type name of field in DB to display: ','s'), go = 1;, end,

while ~isfield(DB,testfield), disp(['NO field called ' testfield]);
    disp(DB), testfield = input('Type field name: ','s');
end
end




