function [OUT] = db_cluster_burt_table(clusters,fnames,DB)
% function [OUT] = db_cluster_burt_table(clusters,field list (cell array of strings),DB)
%
% tor wager
% Prints table of studies and lists whether they found activation in each cluster\
%
% warning: does not count separate CONTRASTS, just STUDIES.
% DB is a structure that contains all info from the database
%
% clusters is a struct in which XYZmm field has list of points
% --should be output of database2clusters, which contains all fields
%
% fnames is a cell array that has names of all fields in clusters to use
% each field should be a list of 1's and 0's in this case, for now.
% enter names in clusters(i).shorttitle to use as region names.
%
% The significance testing done is on the double centered superindicator
% matrix, for an SI matrix with 2 sets of columns: One coding for activation in each cluster, and the other 
% coding for each task type.  Stats are performed on the cross-correlation
% elements between the two sets - i.e., on the associations between
% cluster and task type, controlling for overall number of points in
% clusters and tasks (centering column means), and also controlling for the
% differences between task types in showing more activation
% across all clusters (centering row means for cluster indicators) and
% weighting studies according to the inverse of how many task types they
% involve simultaneously (centering row means for task indicators).

verbose = 0;

disp('Output of db_cluster_burt_table')

[ustudy,suni] = unique(DB.Study);
disp([num2str(length(ustudy)) ' studies in database'])
fprintf(1,'\n')
        
for stud = 1:length(ustudy)

        % get vectors of values in fnames (e.g., 1 or 0) for each study
        % -------------------------------------------------------------------
        for i = 1:length(fnames)
            eval(['tmp1 = DB.' fnames{i} '(suni);'])
            fmatx(stud,i) = tmp1(stud);
        end
            
    
    % get number of peaks in cluster and fill in rmatx (matx of regional peaks)
    % -------------------------------------------------------------------
    for i = 1:length(clusters)
    
        tmp1 = strcmp(ustudy{stud},clusters(i).Study);
        rmatx(stud,i) = sum(tmp1);
     
    end
    
end
        
    
    % Contingency table - header row
    % -------------------------------------------------------------------
    
    print_header(fnames,clusters)
    
    
    % Contingency table - body
    % -------------------------------------------------------------------
    for i = 1:length(ustudy)
        fprintf(1,'%s\t',ustudy{i})
        for j = 1:length(fnames), fprintf(1,'%3.0f\t',fmatx(i,j)), end
        for j = 1:length(clusters), fprintf(1,'%3.0f\t',rmatx(i,j)), end
        fprintf('\n')
    end
  
    
        
    % Compute indicator matrix and Burt table
    % -------------------------------------------------------------------
    
    pmatx = sparse([fmatx rmatx]);
    rmatx = (rmatx > 0);
    bmatx = double([fmatx rmatx]);  % superindicator matrix
    
    % Burt profile plot
    nms = fnames; 
    if isfield(clusters,'shorttitle'), for i = 1:length(clusters), nms = [nms {clusters(i).shorttitle}];, end, end
    for i = 1:size(bmatx,2)-length(nms),nms = [nms {['Reg' num2str(i)]}];, end
    
    percentstudies = burt_profile([bmatx' * bmatx],nms,size(fmatx,2));
    title('Activations by task and brain region')
    
    % create double-centered bmatxc
    tmp1 = double(fmatx); tmp1 = scale(tmp1',1)';
    tmp2 = double(rmatx); tmp2 = scale(tmp2',1)';

    wh1 = find(all(tmp1,2) == 0); 
    wh2 = find(all(tmp2,2) == 0);
    wh = unique([wh1; wh2]);
    
    bmatxc = [tmp1 tmp2];
    bmatxc(wh,:) = [];
    ustudy(wh) = [];
    bmatxc = scale(bmatxc,1);
    
    if ~isempty(wh1)
        disp('Warning: Some rows have the same value for all columns of task (1st set) entries.  Removing.')
    end
    if ~isempty(wh2)
        disp('Warning: Some rows have the same value for all columns of region (2nd set) entries.  Removing.')
    end
    
    burt = scale(bmatxc,1)'*scale(bmatxc,1);    % Burt table of double-centered.  this is what is used in analysis.
    
    
    % Indicator matrix - header row
    % -------------------------------------------------------------------
    fprintf('\n')
    fprintf('Superindicator matrix, row-centered\n')
    print_header(fnames,clusters)
    
    % Indicator matrix - body
    % -------------------------------------------------------------------
    str = repmat('%3.0f\t',1,size(bmatxc,2));
    
    for i = 1:length(ustudy)
        fprintf(1,'%s\t',ustudy{i})
        fprintf(1,str,bmatxc(i,:)), 
        fprintf('\n')
    end
    
    
    %  Burt table - header row
    % -------------------------------------------------------------------
    fprintf('\n')
    fprintf('Burt Table\n')
    print_header(fnames,clusters)
    
    % Burt table - body
    % -------------------------------------------------------------------
    str = repmat('%3.0f\t',1,size(bmatxc,2));
    
    for i = 1:length(fnames)
        fprintf(1,'%s\t',fnames{i})
        fprintf(1,str,burt(i,:)), 
        fprintf('\n')
    end
 
    for i = 1:length(clusters)
        fprintf(1,'%s\t',clusters(i).name)
        fprintf(1,str,burt(i,:)), 
        fprintf('\n')
    end
    
        % Centered Indicator matrix - header row
    % -------------------------------------------------------------------
    %fprintf('\n')
    %fprintf('Double-Centered Superindicator matrix\n')
    %print_header(fnames,clusters)
    
    % Centered Indicator matrix - body
    % -------------------------------------------------------------------
    %str = repmat('%3.6f\t',1,size(bmatx,2));
    
    %for i = 1:length(ustudy)
    %    fprintf(1,'%s\t',ustudy{i})
    %    fprintf(1,str,bmatx(i,:)), 
    %    fprintf('\n')
    %end

    fprintf(1,'\nBootstrapping\n')
    % remove rows so that some studies activating more overall does not influence
    % the correspondence values; remove cols, see above

    OUT = permute_mtx(bmatxc,'ss',1,1000,length(fnames));
    n = length(fnames);
    
    OUT.pmatx = pmatx;
    OUT.descrip_pmatx = 'Sparse superindicator matrix with number of peaks for each study in entries';
    OUT.bmatx = bmatx;
    OUT.descrip_bmatx = 'Superindicator matrix, before centering';
    OUT.n = n;
    OUT.descrip_n = 'First n columns are tasks, rest are regions';
    OUT.ustudy = ustudy;
    OUT.descrip_ustudy = 'Unique studies included';
    OUT.burt = burt;
    OUT.descrip_burt = 'Burt table: Sum of double centered, squared bmatx';
    OUT.percentstudies = percentstudies;
    OUT.descrip_percentstudies = '% of studies activating each region in each task; rows are regions';
    OUT.nms = nms;
    
    % so we have a double-centered matrix; rows and cols are scaled above
    % task x region interaction, double centered
    OUT.MCA = tor_mca(bmatxc,'none',n,nms);
    
    % lines between tasks indicate unequal numbers of studies in diff tasks
    nmdsfig(OUT.MCA.pc,OUT.n,OUT.nms,OUT.p < .05);
    
    % tmp is probability of task, given activity in region (column),
    % without considering prior probabilities of task (e.g., all tasks are
    % equally likely).  
    tmp = OUT.percentstudies' ./ repmat(sum(OUT.percentstudies'),size(OUT.percentstudies,2),1);
    
    
return
    



function print_header(fnames,clusters)

fprintf(1,['\nStudy\t']) 
    for j = 1:length(fnames)    
         fprintf(1,'%s\t',fnames{j}) 
    end

    for j = 1:length(clusters)
        if isfield(clusters,'BAstring'), fprintf(1,'%s\t',clusters(j).BAstring)
        else, fprintf(1,'%s\t',clusters(j).name)
        end
    end
    fprintf('\n\t')
    for j = 1:length(fnames),fprintf(1,'\t'),end
    for j = 1:length(clusters)
        fprintf(1,'%3.0f,%3.0f,%3.0f\t',clusters(j).mm_center(1),clusters(j).mm_center(2),clusters(j).mm_center(3))
    end
    fprintf('\n')
    
return
    