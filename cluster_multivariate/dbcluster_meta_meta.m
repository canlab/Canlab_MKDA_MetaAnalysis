function clnew = dbcluster_meta_meta(varargin)
% clnew = dbcluster_meta_meta(varargin)
%
% Inputs:
%
% clusters structures, one for each meta-analysis
% always followed immediately by DB database structure for the
% meta-analysis
% These structures are the output of database2clusters
%
% cell vector of string myvariable names, one for each meta-analysis
%   - In order of database clusters entered
%   - First cell is myvariable to select, subsequent cells are levels
%
% x, y, and z fields contain coordinates
%
%
% example: 
%dbcluster_meta_meta(EM,EMDB,SW,SWDB,{'valence' 'pos' 'neg'},{'Modality' 'visual'});
% dbcluster_meta_meta(SW,SWDB,PAIN,PAINDB,EM,EMDB,WM,WMDB,{'alltask' 'switching_database'},{'alltask' 'pain_coordinates1'},{'alltask' 'var_emotrev_031602_aav'},{'alltask' 'workingmem'});



% ------------------------------------
% defaults and output files
% ------------------------------------

goglass = 0;

disp('Creating meta_meta_tables.out and meta_meta_chisq.out')

diary meta_meta_tables.out



% ------------------------------------
% concatenate Clusters structures
% ------------------------------------


CL = {}; myvars = {};
clnew = [];

for i = 1:length(varargin)
    if isstruct(varargin{i})
        % clusters structure
        CL{end+1} = varargin{i};
        
    elseif iscell(varargin{i})
        % myvariable names
        
        myvars{end+1} = varargin{i};
        
    end
end

tmp = CL;
CL = tmp(1:2:end);
DBCL = tmp(2:2:end);

% save names of all levels across databases, in correct order.
% Actually shouldn't be necessary to preserve order here...
nms = {};   
colors = {'bo' 'gs' 'r^' 'yv' 'cs' 'mo'};

for i = 1:length(CL)
    
        % Create output structure

        if i == 1
            
            if length(clnew) < 1, 
                clnew.Study = [];
                clnew.conditions = [];
                clnew.rmtx = [];
                clnew.rootn = [];
            end          
        end
        
        % Table of results
        
        clnew.COUNTS(i) = dbcluster_contrast_table(CL{i},DBCL{i},myvars{i});
    
    
        
        % ------------------------------------
        % save meta-analysis data in clusters --> clnew
        % ------------------------------------
        % GET info for this database in this cluster
        
        myvar = myvars{i}(1);                           % name of myvariable to select
        nms = [nms myvars{i}(2:end)];                  % names of levels in order
        
        % list of studies and conditions
        % *** Important step ***
        % Last arg: Turn subject weighting ON
        [fm,rm,ust,taski] = dbcluster2indic(CL{i},DBCL{i},myvar,1); % fm is task list, rm is region mtx
                                                                % taski is
                                                                % task
                                                                % indicator
                                                                % use only
                                                                % to get N
        allfm = []; allrm = []; allust = []; allxyz = []; alltaski = [];
        
        for k = 2:length(myvars{i})
            % find only those studies whose task conditions match desired
            % list
            str = ['wh = find(strcmp(fm, ''' myvars{i}{k} '''));'];
            eval(str);
            
            allfm = [allfm; fm(wh)];        % list of conditions
            allrm = [allrm; rm(wh,:)];      % matrix of region activation indicators
            allust = [allust; ust(wh)];     % unique studies 
            alltaski = [alltaski; taski(wh,:)]; % task indicator with sample size weights
            %allxyz = [allxyz; xyz(wh,:)];   % xyz
        end
        
        % Add to output structure
        
        clnew.Study = [clnew.Study; allust];
        clnew.conditions = [clnew.conditions; allfm];
        clnew.rmtx = [clnew.rmtx; allrm];
        clnew.rootn = [clnew.rootn; sum(alltaski,2)];  %  sample size -- take sqrt below
        % watch out if contrasts belong to more than 1 task!!!  WILL
        % OVERWEIGHT!!!
        
        % ------------------------------------
        % Additional plots of each individual cluster
        % ------------------------------------

        % each region is a different color here
        %regcolors = {'rs' 'rs' 'yd' 'yd' 'bv' 'bv' 'gs' 'gs'};
        %dbcluster_glass_plot(CL{i},regcolors);
        %saveas(gcf,CL{i}(1).alltask{1},'fig')
        %saveas(gcf,CL{i}(1).alltask{1},'tif')
        
        % if we have multiple conditions within, plot these
        
        %dbcluster_glass_plot(CL{i},colors,myvar{1},nms);
        %saveas(gcf,[CL{i}(1).alltask{1} '_levels3'],'fig')
        %saveas(gcf,[CL{i}(1).alltask{1} '_levels3'],'tif')
        %f = gcf; figure(f-1)
        %saveas(gcf,[CL{i}(1).alltask{1} '_levelsR'],'fig')
        %saveas(gcf,[CL{i}(1).alltask{1} '_levelsR'],'tif')
        
        clnew.CL{i} = CL{i};
        
end % database

clnew.rootn = clnew.rootn .^ .5;    % square root

diary off

% --------------------------------------------------------------
% Glass plot
% --------------------------------------------------------------
i = 1;

if goglass,
dbcluster_glass_plot(CL{i},colors(i),myvar{1},myvars{i}(2:end));
for i = 2:length(CL)
    dbcluster_glass_plot(CL{i},colors(i),myvar{1},myvars{i}(2:end),'add'); drawnow
end
saveas(gcf,'meta_meta_glass3','fig')
saveas(gcf,'meta_meta_glass3','tif')
f = figure(gcf);
figure(f-1);
saveas(gcf,'meta_meta_glassR','fig')
saveas(gcf,'meta_meta_glassR','tif')
end % if go

% --------------------------------------------------------------
% Bar plot
% --------------------------------------------------------------

studycount = []; conditioncount = []; levels = {};
for i = 1:length(clnew.COUNTS),
    studycount = [studycount; clnew.COUNTS(i).cond_by_cl_counts];
    conditioncount = [conditioncount; clnew.COUNTS(i).conditioncounts];
    levels = [levels clnew.COUNTS(i).levels{:}];
end
out = perc_barplot(studycount,conditioncount,levels,'meta_meta',clnew.COUNTS(1).clusternames);
barfig = gcf;

domontage = 0;

if domontage
    
if ~isempty(out.AvsB) & ~isempty(out.BvsA) 
    montage_clusters([],CL{1}(out.AvsB),CL{1}(out.BvsA),colors);
elseif ~isempty(out.AvsB) 
     montage_clusters([],CL{1}(out.AvsB),colors);
elseif ~isempty(out.AvsB) 
     montage_clusters([],CL{1}(out.BvsA),colors(2));
elseif ~isempty(out.Greatest)
    % more than 2 conditions
    u = unique(out.Greatest)';
    gcol = colors(u);
    str = ['montage_clusters([],'];
    for i = u    % find significant regions
        str = [str 'CL{1}(find(out.Greatest==' num2str(i) ')),'];
    end
    str = [str 'gcol);'];
    disp(str)
    eval(str)
    
else
    
    disp('No significant regions');
end

saveas(gcf,['montage_' clnew(1).COUNTS(1).factornames{1}],'fig');
saveas(gcf,['montage_' clnew(1).COUNTS(1).factornames{1}],'tif');

end % if domontage



% On new clusters, get indicator matrix, compiled across all clusters
% This is relevant when more than 1 database structure is entered
% clnew.indic is indicators over all tasks in all databases, + regions
% clnew.task_indicator is tasks only
% clnew.rmatx is peak locations
conds = clnew.conditions;
for i = 1:length(conds), conds{i}(find(conds{i}=='_')) = ' ';, end
[clnew.indic,clnew.nms] = string2indicator(conds,nms); 

clnew.rootn(find(isnan(clnew.rootn))) = nanmean(clnew.rootn);
clnew.indic = clnew.indic .* repmat(clnew.rootn,1,size(clnew.indic,2));

% unnecessary now
% replace NaNs in case of missing subjects
%for i = 1:size(clnew.indic,2), 
%    wh = find(isnan(clnew.indic(:,i)));
%    clnew.indic(wh,i) = nanmean(clnew.indic(:,i));
%end
clnew.task_indicator = clnew.indic(:,1:size(unique(clnew.conditions),1));   % first k cols, k is number of factors

% get names for regions

if isfield(CL{1},'shorttitle'), for i = 1:length(CL{1}), nms = [nms {CL{1}(i).shorttitle}];, end, end
sz = size(clnew.rmtx,2) + size(clnew.indic,2) -length(nms);
for i = 1:sz,nms = [nms {['Reg' num2str(i)]}];, end

for i = 1:length(nms), nms{i}(nms{i} == '_') = ' ';,end
clnew.nms = nms;
clnew.n = size(clnew.indic,2);
clnew.bmatx = [clnew.indic clnew.rmtx];




% --------------------------------------------------------------
% Add p(T|A) predictors to barplot
% --------------------------------------------------------------

figure(barfig)

if length(clnew(1).COUNTS) == 1 % only 1 database entered
    [ind,col] = pt_given_a_plot(clnew.COUNTS.effpt_given_a,clnew.COUNTS.Ppt_given_a,levels,'PTplot_',clnew.COUNTS(1).clusternames);
    pt_given_a = clnew.COUNTS.pt_given_a;
else
    for i = 1:size(clnew.rmtx,2)     % for each region; weighting ON
        [pt_given_a{i},eff{i},z{i},p{i},stat] = p_task_given_activation(clnew.task_indicator,clnew.rmtx(:,i),1); 
        pt_given_a{i} = pt_given_a{i} ./ sum(pt_given_a{i});
    end
    [ind,col] = pt_given_a_plot(eff,p,levels,'PTplot_',clnew.COUNTS(1).clusternames);
end

saveas(gcf,['PTplot_' clnew(1).COUNTS(1).factornames{1}],'fig');
saveas(gcf,['PTplot_' clnew(1).COUNTS(1).factornames{1}],'tif');

% likelihood ratio, pt | a / pt
% to remove effects of task frequency -- flat priors on task
% normalize across tasks 
v = cat(1,pt_given_a{:});
clnew.pt = sum(clnew.task_indicator) ./ sum(clnew.task_indicator(:));
v = v ./ repmat(clnew.pt,size(v,1),1);
v = v ./ repmat(sum(v,2),1,size(v,2));

stackedbar2(v,[],nms)
set(gca,'XLim',[0 1])
ylabel('Brain region'),title('Predictions of task given brain activity')
xlabel('Probability') 
set(gcf,'Position',[1363         621        1053         490])
saveas(gcf,[clnew(1).COUNTS(1).factornames{1} '_predict_single'],'fig')
saveas(gcf,[clnew(1).COUNTS(1).factornames{1} '_predict_single'],'tif')



% classify based on effects
% conditional on finding some activation somewhere

v = cat(1,eff{:});

% studies that fit well vs. poorly: v * rmatx'
% but watch out for N weighting

w = clnew.rmtx * v;                     % weights on each task for each study

wh = find(all(w == 0,2));
fprintf(1,'%3.0f contrasts activated no regions: no info to classify',length(wh));
w(wh,:) = [];

w = w - repmat(max(w')',1,size(w,2));   % zeros where max
[R,C] = find(w == 0);                    % C is classes
% worry about ties later

% true
[Rtrue,Ctrue] = find(clnew.task_indicator);
Rtrue(wh) = [];
Ctrue(wh) = [];
[m,dp,corr,far] = confusion_matrix(Ctrue,C)
% worry about non-unique task classes later


% try another classification
rmtx = clnew.rmtx; rmtx(wh,:) = [];
rmtx = real(rmtx>0);

for i = 1:size(v,2)     % for each task class
    tmp = nanmean(rmtx(Ctrue == i,:));
    truemean(i,:) = tmp;
end
%truemean = truemean ./ repmat(sum(truemean,2),1,size(truemean,2));
%rmtx = real(rmtx>0) ./ repmat(sum(rmtx,2),1,size(rmtx,2));
for i = 1:size(rmtx,2)
    d = repmat(rmtx(i,:),size(truemean,1),1);
    d = sum((d - truemean).^2,2);
    C(i) = find(d == min(d));
end
[m,dp,corr,far] = confusion_matrix(Ctrue,C)


% not working so hot
%tor_fig;
%ind = ind';
%for i = 1:length(ind), 
%    subplot(ceil(length(ind)/3),3,i); 
%    hold on; montage_clusters_maxslice([],CL{1}(i),col{i});,
%end

%x = left+3*(.8/21); y = bottom + .65;
%h = axes('position',[x y .95/21 .15]);
%montage_clusters_maxslice([],CL{1}(3),{'r'});

clnew.chi2 = out;


% clnew is structure
% has one row per contrast
% COUNTS(i) has info about each separte database
% indic has indicators for conditions (task type)
% rmatx has data about activated regions (counts per region, by contrast
% within study)
% bmatx has both




doclassify = 1;
if doclassify

% and classify them
diary meta_meta_tables.out

[clnew.prob clnew.perc] = indic2classify1(clnew.bmatx,clnew.n,nms);
nms
clnew.OUT = indic2classify2(clnew.bmatx,clnew.n,nms);

% show which studies were misclassified
fprintf(1,'\nKNN Misclassified Contrasts\n');
fprintf(1,'\tClassed as\t');
for i=1:length(nms), fprintf(1,'%s\t',nms{i});, end, fprintf(1,'\n');
wh = find(clnew.OUT.misclass);
for i = 1:length(wh),
    fprintf(1,'%s\t',clnew.Study{wh(i)})
    fprintf(1,'%s\t',nms{clnew.OUT.class(wh(i))})
    fprintf(1,'%3.0f\t',clnew.bmatx(wh(i),:));
    fprintf(1,'\n');
end


diary off

end


dbclusters = clnew;
str = ['save dbcluster_' clnew(1).COUNTS.factornames{1} ' dbclusters'];
eval(str);

return


        
    
        
        
        
        %indic = zeros(size(xyz,1),length(myvars{i})-1);
        %for k = 2:length(myvars{i})
            
        %    str = ['tmp = strcmp(CL{i}(j).' myvars{i} ' == ' myvars{k} ');'];
        %    eval(str);
        %    indic(:,k-1) = tmp;
        %    nms{k-1} = [myvar '_' myvars{k}];
            
        end
        
        