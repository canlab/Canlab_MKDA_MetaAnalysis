function outtable = cluster_chi2(clusters,fnames,DB)
% function outtable = cluster_chi2(clusters,field list (cell array of strings),DB)
%
% tor wager
% counts studies in each cluster
% warning: does not count separate CONTRASTS, just STUDIES.
% DB is a structure that contains all info from the database
%
% clusters is a struct in which XYZmm field has list of points
% --should be output of database2clusters, which contains all fields
%
% fnames is a cell array that has names of all fields in clusters to use
%
% assume field list variables code for two tasks, A and B.
% out of the set of studies that reported peaks in a cluster (brain
% region), what proportion of studies of A produced activations in the
% region?  What proportion of studies of B did?  Are the two proportions
% different?  The chi2 measure generalizes this across N conditions.
%
% Note: this function does not account for the total number of activations
% produced by studies of a particular type.  e.g., If one task is a better
% 'activator' overall, then it will produce significant results in all
% clusters tested.
%
% This function also does not test the specificity of activations in a
% brain region

% scount    : studycount
%   columns are Yesses, Nos.  Rows are conditions
%   within one cluster (region)
%
% pyes      : prob of yes
%   expected proportion of studies activating a region compared to total
%   studies.  expectation is based on average across conditions.
%   = sum(scount(:,1)) ./ sum(stotals)
%
% totals
%   sums of yesses and nos (sum(scount,2))
% 

verbose = 0;

disp('Output of cluster_chi2 function')

[tmp,suni] = unique(DB.Study);
disp([num2str(length(tmp)) ' studies in database'])
fprintf(1,'\n')
        
if ~isfield(clusters,'XYZall'), 
    for i = 1:length(clusters)
        clusters(i).XYZall = [clusters(i).x clusters(i).y clusters(i).z];, 
    end
end

docellinput = 0;
if length(fnames) == 1 & iscell(DB.(fnames{1})(1))
    disp('Input field is names; converting to indicators');
    docellinput = 1;
    
    [fmatx,rmatx,ustudyout,taskindic,tasknms,allsuni] = dbcluster2indic(clusters,DB,fnames);
    whomit = sum(taskindic,2); whomit = find(whomit==0);  % find NaN entries in neither task and remove
    taskindic(whomit,:) = [];
    
    bnull = sum(taskindic) ./ size(taskindic,1);
    stotals = sum(taskindic);
    
    for j = 1:length(tasknms)
        disp([tasknms{j} ': ' num2str(stotals(j)) ' studies, ' num2str(100*bnull(j)) '% of studies'])   
    end
    

else
    % Count OVERALL did vs. did not involve this condition
    % -------------------------------------------------------------------
    for j = 1:length(fnames)

            str = ['tmp12 = DB.' fnames{j} ';'];
            eval(str)



            tmp12 = tmp12(suni);
            bnull(j) = sum(tmp12 == 1) ./ length(tmp12);

            stotals(j) = sum(tmp12 == 1);
            disp([fnames{j} ': ' num2str(stotals(j)) ' studies, ' num2str(100*bnull(j)) '% of studies'])
    end
    
end



if docellinput
    fnames = repmat(fnames,1,length(tasknms));
end
    

for i = 1:length(clusters)
    
    numstudies = length(unique(clusters(i).Study));
    numpts = size(clusters(i).XYZall,1);
    disp(['Cluster ' num2str(i) ': ' num2str(numpts) ' points in ' num2str(numstudies) ' studies.'])
    
   
    for j = 1:length(fnames)
        
        % count Yesses for this condition
        % -------------------------------------------------------------------
        if docellinput
            stud = clusters(i).Study(strcmp(clusters(i).(fnames{1}),tasknms{j}));

        else
            % if field contains indicators (1/0)
            str = ['stud = clusters(i).Study(find(clusters(i).' fnames{j} '));'];
            eval(str)
        end
        
        [allstud,uni] = unique(clusters(i).Study);
        [ustud] = unique(stud);
        scount(j,1) = length(ustud);
            
        if verbose
            fprintf(1,'\tStudies involving %s\n',fnames{j})
        
            for k = 1:length(ustud)
            
                pks(k) = sum(strcmp(stud,ustud{k}));
                fprintf(1,'\t\t%s %3.0f pk\n',ustud{k},pks(k))
            
            end
            
        end
        
        % count Nos for this condition
        % -------------------------------------------------------------------
        
        if docellinput
            stud = clusters(i).Study(~strcmp(clusters(i).(fnames{1}),tasknms{j}));
        else
            str = ['stud = clusters(i).Study(find(~clusters(i).' fnames{j} '));'];
            eval(str)
        end
        
        ustud = unique(stud);
        scount(j,2) = length(ustud);
                 
        if verbose
            fprintf(1,'\tStudies NOT involving %s\n',fnames{j})
        
            for k = 1:length(ustud)
            
                fprintf(1,'\t\t%s %3.0f pk\n',ustud{k},sum(strcmp(stud,ustud{k})))
            
            end
        
            fprintf(1,'\n')
        end
               
        % save overall criterion column
        % -------------------------------------------------------------------
        str = ['tmp11 = clusters(i).' fnames{j} ';'];
        eval(str)
        
        if docellinput
            x{i}(:,j) = double(strcmp(tasknms{j},tmp11(uni)));
        else
            x{i}(:,j) = tmp11(uni);
        end
        
    end
    
    % Chi2 and binomial tables
    % -------------------------------------------------------------------
    yesses = scount(:,1)';
    totals = sum(scount,2)';
    [chi2,colnms] = computechi2(yesses,stotals);     %totals); 
    fprintf(1,'\n\tOmnibus chi-square test for diffs among conditions in p of study reporting peak\n\t')

    fprintf(1,'Chi2\tdf\tp\tCramV\tSig?\tWarn\tPyes\t\n'), 

    fprintf(1,'\n\t'), fprintf(1,'%3.3f\t',chi2)
    fprintf(1,'\n')
    
    if i == 1, outtable = sprintf('Chi2\tdf\tp\tCramV\tSig?\tWarn\tPyes\t\n');, ,end
    outtable = str2mat(outtable,sprintf('\t',chi2));
    outtable = str2mat(outtable,sprintf('%3.3f\t',chi2));
     
    b = binocdf(yesses,totals(1),bnull);          % chance prob of same or fewer than observed hits 
    b2 = (1-b) + binopdf(yesses,totals(1),bnull); % chance prob of same or more than observed hits
    
    %b = binopdf(yesses,totals(1),bnull);
    %b = min(b,1-b);
    fprintf(1,'\n\tP-values for binomial test within each condition\n')
    fprintf(1,'\t(test that proportions of studies in this condition is reflected in activations in this region.\n')
    fprintf(1,'\te.g., if 16%% of studies involved condition A, then 16%% of activations in R involve condition A)\n')
    
    for j = 1:length(fnames), fprintf(1,'\t\t%s ',fnames{j}),end
    fprintf(1,'\nProb. that fewer than expected activations are observed\t'), fprintf(1,'%3.3f\t',b)
    fprintf(1,'\nProb. that more than expected activations are observed\t'), fprintf(1,'%3.3f\t',b2)
    fprintf(1,'\n')
    
    % Contingency table
    % -------------------------------------------------------------------
    fprintf(1,['\n\tStudy\t']) 
    for j = 1:length(fnames)    
        if docellinput
            fprintf(1,'%s\t',tasknms{j}) 
        else
            fprintf(1,'%s\t',fnames{j}) 
        end
    end
    fprintf('\n\t')
    for j = 1:length(allstud)
        fprintf(1,['%s\t' repmat('%3.0f\t',1,length(fnames)) '\t\n\t'],allstud{j},x{i}(j,:))
    end
    fprintf(1,['\n\tHit rate\t']) 
    for j = 1:length(fnames)    
         fprintf(1,'%3.2f\t',yesses(j) ./ stotals(j)) 
    end
    fprintf('\n\t')
    
end


return
    