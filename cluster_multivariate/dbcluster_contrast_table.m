function OUT = db_cluster_table(clusters,DB,varargin)
% function OUT = db_cluster_table(clusters,DB,varargin)
%
% tor wager
% counts studies and contrasts in each cluster
%
% clusters is a struct in which XYZmm field has list of points
% --should be output of database2clusters, which contains all fields
% DB is database output of database2clusters
%
% Prints one row per coordinate in any cluster, along with study names,
% cluster ID, and values in fields you specify in fnames
%
% optional arguments:
% fnames is a cell array that has names of all fields in clusters to use
% first cell is which factor, following cells are levels to include
% enter one argument per factor:
%
% example:
% OUT = dbcluster_contrast_table(PAIN,PAINDB,{'Right_vs_Left' 'Right' 'Left'})
% perc_barplot(OUT.cond_by_cl_counts,OUT.conditioncounts,OUT.levels{1},'test',OUT.clusternames)
%
% see also db_cluster_table.m 

% do taskindicators
dotaskindic = 0;


% enter input arguments

for i = 1:length(varargin), 
    fnames{i} = varargin{i}{1};
    if length(varargin{i}) > 1
        levels{i} = varargin{i}(2:end);
    else
        eval(['levels{i} = unique(DB. ' varargin{i}{1} ');'])
    end
end


% get indicator for all levels
% Sample Size weighting ON
[fmatx,rmatx,ustudyout,taskindic,tasknms,allsuni] = dbcluster2indic(clusters,DB,fnames,1);

fprintf(1,'%3.0f Unique contrasts\n',length(ustudyout))

% select only levels of interest
newti = []; newnames = {};

for i = 1:length(fnames)
    for j = 1:length(levels{i})
        wh = find(strcmp(tasknms,levels{i}{j}));
        
        newti = [newti taskindic(:,wh)];
        newnames = [newnames tasknms(wh)];
    end
end
           

% Display header

fprintf(1,'Study\t')
if dotaskindic
for i = 1:length(newnames)
    fprintf(1,'%s\t',newnames{i})
end
end

    for i = 1:size(fnames,2)
        fprintf(1,['%s\t'],fnames{i})
    end
    
    for i = 1:size(rmatx,2)
        if ~isfield(clusters,'shorttitle')
            fprintf(1,['Clust%3.0f\t'],num2str(i))
        else
            fprintf(1,['%s\t'],clusters(i).shorttitle)
        end
    end
    
fprintf(1,'\n')


% Display table body
[yr,ustudyout] = getYear(ustudyout);

for i = 1:length(ustudyout)
    
    
    % print study name, task values, and number of peaks in each cluster
    fprintf(1,'%s\t',[ustudyout{i} ' ' yr{i}]);
    
    if dotaskindic
    for j = 1:size(newti,2)
        fprintf(1,['%3.0f\t'],newti(i,j))
    end
    end

    for j = 1:size(fmatx,2)
        fprintf(1,['%s\t'],fmatx{i,j})
    end
    
    for j = 1:size(rmatx,2)
        fprintf(1,['%3.0f\t'],rmatx(i,j))
    end
    
    fprintf(1,'\n')
    
end

% Summary and percentages for columns (brain regions)

totals = sum(newti > 0);            % task totals
clustertotals = sum(rmatx > 0);     % overall activation in each region

if ~dotaskindic         % repeat headers to avoid confusion
    fprintf(1,'\n\t\t')               % task indicator totals
    for i = 1:length(newnames)
        fprintf(1,'%s\t',newnames{i})
    end
    fprintf(1,'\n\t')               % cluster count totals
    fprintf(1,['Total\t' repmat('%3.0f\t',1,size(totals,2))],totals);
    fprintf(1,'\n\t\t')
    for i = 1:size(rmatx,2)
        if ~isfield(clusters,'shorttitle')
            fprintf(1,['Clust%3.0f\t'],num2str(i))
        else
            fprintf(1,['%s\t'],clusters(i).shorttitle)
        end
    end
    fprintf(1,'\n\t\t')  
    fprintf(1,[repmat('%3.0f\t',1,size(clustertotals,2))],clustertotals);
else        

    fprintf(1,'\n')
    fprintf(1,['Total\t' repmat('%3.0f\t',1,size(totals,2))],totals);
    for i = 1:size(fnames,2), fprintf(1,'\t'),end

    fprintf(1,[repmat('%3.0f\t',1,size(clustertotals,2))],clustertotals);
end

% percentages for columns (brain regions)

perc = 100* sum(newti > 0) ./ length(ustudyout);
clusterperc = 100 * sum(rmatx > 0) ./ length(ustudyout);

fprintf(1,'\n')
fprintf(1,['Percentage\t' repmat('%3.0f\t',1,size(perc,2))],perc);
for i = 1:size(fnames,2), fprintf(1,'\t'),end

if ~dotaskindic         % line break if no task indicator
    fprintf(1,'\n')
end

fprintf(1,[repmat('%3.2f\t',1,size(clusterperc,2))],clusterperc);

fprintf(1,'\n')
fprintf(1,'\n')


% replace NaN task values with 0
wh = find(isnan(newti)); newti(wh) = 0;


% Summary and percentages by condition


for j = 1:size(newti,2) % newti is task indicator
    
    for k = 1:size(rmatx,2)
        
        clcondtotal(j,k) = sum(real(rmatx(:,k) > 0) & newti(:,j));
        clcondperc(j,k) = 100 * clcondtotal(j,k) / totals(j);   % % of studies in condition
        
    end
    
end


% Prob of task given activation in each region 


% p(t|a) = p(a|t)*p(t) / p(a)

% with flat priors -- old way
%pa = clustertotals ./ length(ustudyout);    % overall p(a) for each region
%pt = ones(1,size(newti,2)) ./ size(newti,2);    % overall flat priors on task: 1 / number of tasks
%pt = repmat(pt',1,size(clustertotals,2));       % replicate for each region; all rows of pt should be the same 
%pta = clcondperc ./ 100;                        % p(a|t) estimate
%for j = 1:size(newti,2)     % for each task condition
%    pt_given_a(j,:) = pta(j,:) .* pt(j,:) ./ pa;   
%end

for i = 1:size(rmatx,2)     % for each region; weighting ON
    [pt_given_a{i},eff{i},z{i},p{i},stat] = p_task_given_activation(newti,rmatx(:,i),1); 
end




        % ------------------------------------
        % print tables
        % ------------------------------------
        
        
fprintf(1,'Totals by condition and cluster\n')

% header

fprintf(1,'\t')
    for i = 1:size(rmatx,2)
        if ~isfield(clusters,'shorttitle')
            fprintf(1,['Clust%3.0f\t'],num2str(i))
        else
            fprintf(1,['%s\t'],clusters(i).shorttitle)
        end
    end
    fprintf(1,'\n')
    
for j = 1:size(newti,2)
    
    fprintf(1,'%s\t',newnames{j})
    
    str = [repmat('%3.0f\t',1,size(clcondtotal,2))];
    fprintf(1,str,clcondtotal(j,:))
    fprintf(1,'\n')
    
end


fprintf(1,'Percentages by condition and cluster\n')
for j = 1:size(newti,2)
    
    fprintf(1,'%s\t',newnames{j})
    
    str = [repmat('%3.0f\t',1,size(clcondperc,2))];
    fprintf(1,str,clcondperc(j,:))
    fprintf(1,'\n')
    
end


fprintf(1,'Probability of task given activation in each region\n')
fprintf(1,'Weighting: %s\n',stat.weighted)
       
for i = 1:size(rmatx,2)     % regions
        if ~isfield(clusters,'shorttitle')
            fprintf(1,['Clust%3.0f\t'],num2str(i))
        else
            fprintf(1,['%s\t'],clusters(i).shorttitle)
        end

        % header for task conditions
    fprintf(1,'\t')
    for j = 1:size(newti,2)     % task indicators of interest   
        fprintf(1,'%s\t',newnames{j})
    end
    fprintf(1,'\n')

    fprintf(1,'Pt|A - Pt\t');
    % effect pt_given_a - pt
    str = [repmat('\t%3.2f\t',1,length(eff{i}))];
    fprintf(1,str,eff{i})
    fprintf(1,'\n')
    
    fprintf(1,'Z\t');
    % z-score of effect pt_given_a - pt
    fprintf(1,str,z{i})
    fprintf(1,'\n')
    
    fprintf(1,'p\t');
    % z-score of effect pt_given_a - pt
    fprintf(1,str,p{i})
    fprintf(1,'\n')
    
    fprintf(1,'\n')
end
    


        % ------------------------------------
        % save output in OUT structure
        % ------------------------------------
        
OUT.cond_by_cl_counts = clcondtotal;
OUT.clustercounts = clustertotals;
OUT.conditioncounts = totals;
OUT.ustudyout = ustudyout;
OUT.fmatx = fmatx;
OUT.rmatx = rmatx;
OUT.factornames = fnames;
OUT.levels = levels;
OUT.pt_given_a = pt_given_a;
OUT.effpt_given_a = eff;
OUT.Zpt_given_a = z;
OUT.Ppt_given_a = p;

for i = 1:length(clusters),
    if isfield(clusters,'shorttitle')
        OUT.clusternames{i} = clusters(i).shorttitle;,
    else
        OUT.clusternames{i} = ['CL' num2str(i)];
    end
end

%perc_barplot(OUT.cond_by_cl_counts,OUT.conditioncounts,OUT.levels{1},'test',OUT.clusternames)
  
return
  
  
  

function [b,a] = getYear(ustudy)
% b is year, a is study
clear a, clear b
for i = 1:length(ustudy)
    a{i} = ustudy{i}(1:end-2);
    b{i} = ustudy{i}(end-2:end);
    
    if strcmp(b{i}(end-1:end),'00') | strcmp(b{i}(end-1:end),'01') | strcmp(b{i}(end-1:end),'02') | ...
            strcmp(b{i}(end-2:end-1),'00') | strcmp(b{i}(end-2:end-1),'01') | strcmp(b{i}(end-2:end-1),'02')
        c = '20';
    else
        c = '19';
    end
    
    if ~(strcmp(b{i}(end),'a') | strcmp(b{i}(end),'b') | strcmp(b{i}(end),'c'))
        b{i} = b{i}(end-1:end);
    end
    b{i} = [c b{i}];
end



return

