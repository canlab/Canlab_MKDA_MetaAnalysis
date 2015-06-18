function [fmatx,rmatx,ustudyout,taskindic,tasknms,allsuni] = dbcluster2indic(clusters,DB,fnames,varargin)
%[fmatx,rmatx,ustudyout,taskindic,tasknms,indices (allsuni)] = dbcluster2indic(clusters,DB,fnames,[Sample Size])
%
%[level ID, point count, study ID, indicator of study (contrast) x levels,unique level names] = ...
%
% makes an indicator matrix, given clusters structure and database
% structure (containing all database info)
% fnames is a cell array that has names of all fields in clusters to use
% 
% INPUTS:
%
% optional: sample size flag
%
% OUTPUTS:
%
% fmatx matrix of cell strings, contrasts w/i study=rows,
%   conditions/factors=columns, entries are levels of conditions
% rmatx matrix peaks in contrasts=rows, clusters=columns
% ustudyout is list of study names for unique contrasts, with 1 study name
%   per unique contrast within study
% taskindic : indicator matrix of task (contrast) space, 1/0 if the contrast does or
%   does not involve each level of each factor.
% tasknms : names for each column of taskindic in cell vector of strings
%
% allsuni   : index of which rows in DB are unique contrasts, so that
%   ustudyout = DB.Study(allsuni);
%
% Example:
% [fm,rm,ust] = dbcluster2indic(EM,EMDB,{'valence'})
%
% 
% COUNTS CONTRASTS!  this is good.

% ------------------------------------
% initialization and checks
% ------------------------------------

fmatx = []; rmatx = [];ustudyout = []; taskindic = [];tasknms = [];allsuni =[];

% change fnames (var names) to cell if it's not
if ~iscell(fnames), tmp=fnames; fnames={}; for i =1:size(tmp,1),fnames{i}=tmp(i,:);,end,end

% check to see if this is a clusters structure with DB info
if ~isfield(clusters,'Study')
    if isfield(clusters,'study')
        for i = 1:length(clusters), clusters(i).Study = clusters(i).study;,end
    else
        error('No Study variable.  Does cl have meta-database info?  Use database2clusters first.');
    end
end

[ustudy,suni] = unique(DB.Study);

% preserve order
[suni,wh]=sort(suni); ustudy=ustudy(wh); 

disp([num2str(length(ustudy)) ' studies in database'])
fmatx = {};
allsuni = [];       % all rows extracted as unique contrasts
taskindic = [];     % keep track of task levels for each unique contrast in each study


% ------------------------------------
% Get levels of factor
% ------------------------------------

% get all levels of all factors
for i = 1:length(fnames)
    eval(['[levels{i},x1] = unique(DB.' fnames{i} ');'])
    
    % preserve order
    [x2,whx]=sort(x1); levels{i}=levels{i}(whx);

    if ~iscell(levels{i}), 
        levels{i} = levels{i}(~isnan(levels{i})); tmp12 = {};
        for ii = 1:length(levels{i}), tmp12{ii,1} = num2str(levels{i}(ii));,end, levels{i} = tmp12;,
    else
        % remove NaNs, cause they cause problems
        whom = [];
        for ii = 1:length(levels{i}),
            if isnan(str2num(levels{i}{ii})), whom(ii) = 1;, else, whom(ii) = 0;,end
        end
        levels{i}(find(whom)) = [];
    end
    
end

for stud = 1:length(ustudy)

    thisstudy = find(strcmp(DB.Study,ustudy{stud}));    % all points in this study
    tmpindic = [];                                      % for counting unique contrasts
    e = size(fmatx,1);                                  % last elem in fmatx
    taskindictmp = [];                                     % for storing indicators for task levels
    tasknms = [];                                       % names of all levels of all factors
    
        % get vectors of values in fnames (e.g., 1 or 0) for each study
        % in other words, get an indicator matrix of points (rows) x task conditions for the study
        % get unique combinations of levels of fnames - these are contrasts
        % in the study
        % -------------------------------------------------------------------
        for i = 1:length(fnames)    % for each variable analyzed
            eval(['tmp1 = DB.' fnames{i} '(thisstudy);'])   % tmp1 is condition of each point
            
            [tmpindic] = [tmpindic string2indicator(tmp1)]; % has indicators for all levels of all factors
                            % will return a matrix if more than one unique
                            % value in tmp1; doesn't necessarily respect
                            % levels

            % keep track of overall levels of each factor in separate
            % indicator
            [ti,tasknmstmp] = string2indicator(tmp1,levels{i}); % ti: points in study x levels indicator
                                                             % tasknmstmp:
                                                             % names of
                                                             % levels
                                                              
                                                               
            taskindictmp = [taskindictmp ti];       % taskindictmp = ALL factors, indicator of all levels
            tasknms = [tasknms tasknmstmp];         % names of all levels of all factors
        end
        
        % sunistudy is vector for unique contrasts within THIS study
        % containing index numbers in overall DB for unique contrasts in
        % this study
        [uu,sunistudy] = unique(tmpindic,'rows');           % unique contrasts
        
        
        taskindictmp = taskindictmp(sunistudy,:);           %  task combos for unique contrasts in this study
        
        sunistudy = thisstudy(1) + sunistudy - 1;   % convert to overall index in DB
        len = length(sunistudy);
        allsuni = [allsuni; sunistudy];             % unique contrasts concatenated for all studies
        
        % -----------------------------------------------------------------
        % add Subject Count (sample size) to task indicator here, if
        % requested
        % -----------------------------------------------------------------

            if length(varargin) > 0
                subs = DB.Subjects(sunistudy);
                if size(subs,1) > size(subs,2), subs = subs';, end  % make it a row vector
                % if empty, use mean subs across studies
                if iscell(subs)
                  if isempty(subs{1})
                    subs = DB.Subjects;
                    % convert to cell if it's not; assume string if cell
                    clear subs2
                    if iscell(subs), 
                        for si = 1:length(subs), 
                            if ~isempty(subs{si})
                                subs2(si)=str2num(subs{si});,
                            else
                                subs2(si) = NaN;
                            end
                        end; subs = subs2;, 
                    end
                    
                    subs = nanmean(subs);
                    subs = repmat(subs,1,size(taskindictmp,1));
                  end
                end
                    
                % convert to cell if it's not; assume string if cell
                clear subs2
                if iscell(subs), 
                    for si = 1:length(subs), subs2(si)=str2num(subs{si});,end; subs = subs2;, 
                end
            
                subs = repmat(subs',1,size(taskindictmp,2));
                % save sample sizes in task indic
                taskindic = [taskindic; taskindictmp.* subs];  
                
            else
        
                taskindic = [taskindic; taskindictmp];  % save indicator of task levels for these contrasts
            end
            
        
        tmpindic = [];

        % -----------------------------------------------------------------
        % get indicator with rows = unique contrasts, for building rows of
        % output indicator 
        % -----------------------------------------------------------------
        
        for i = 1:length(fnames)    % loop through levels * factors of interests
            eval(['tmp1 = DB.' fnames{i} '(sunistudy);'])
            % tmp1 is array of values for this particular task type for
            % this study
            % make cell array if it's not
            if ~iscell(tmp1), clear tmp12, for ii = 1:length(tmp1), tmp12{ii,1} = num2str(tmp1(ii));,end, tmp1 = tmp12;,end
            
            fmatx(e+1:e+len,i) = tmp1;  % matrix of level names (rows) for each factor (cols) for unique contrasts
            
            %ustudyout(e+1:e+len) = ustudy(stud);   % done at end
            
            [tmp2,tmpnms{i}] = string2indicator(fmatx(e+1:e+len,i));
            
          
            % -----------------------------------------------------------------
            % add indicator for these unique contrasts for this study
            % (tmp2) to all indicators in this study
            % -----------------------------------------------------------------           
            % tmpindic is indicator for unique contrasts w/i this study???
            tmpindic = [tmpindic tmp2];  % save for later matching to clusters
            
        end
        
        studyindic = [];
        % get indicator with rows = points in this study, for matching rows
        % to cluster values and getting number of peaks of each contrast
        % type in each cluster
        for i = 1:length(fnames)
            eval(['tmp1 = DB.' fnames{i} '(thisstudy);'])
            studyindic = [studyindic string2indicator(tmp1)];  % save for later matching to clusters
        end
        


    % get number of peaks in cluster and fill in rmatx (matx of regional peaks)
    % -------------------------------------------------------------------
    for i = 1:length(clusters)
    
        thisstudy = find(strcmp(ustudy{stud},clusters(i).Study));
        tmpindic2 = [];
        
        if isempty(thisstudy)
            rmatx(e+1:e+len,i) = 0;         % matrix peaks in contrasts=rows, clusters=columns
        else
            % find indicator matrix for this cluster, with columns that match earlier indicator matrix for this study 
            % first loop through factors
            for j = 1:length(fnames)
                eval(['tmp1 = clusters(i).' fnames{j} '(thisstudy);'])
            
                if ~iscell(tmp1), clear tmp12, for ii = 1:length(tmp1), tmp12{ii,1} = num2str(tmp1(ii));,end, tmp1 = tmp12;,end
                
                for k = 1:length(tmpnms{j}) % for each level of this factor
                    tmpindic2 = [tmpindic2 strcmp(tmp1,tmpnms{j}{k})];
                end
                 
                % this won't work because it doesn't always preserve factor
                % levels from the original tmpindic, the indicator for the
                % study across all clusters
                %[tmpindic2] = [tmpindic2 string2indicator(tmp1)];
            end
        
            % go through each contrast in study, and find how many
            % in-cluster points (rows of indicator) match the indicator row
            % for that contrast
            
            for k = 1:size(tmpindic,1)  % for each contrast
                
                % tmpindic is the indicator for the unique contrast
                % tmpindic2 is for each point
                
                % build indicator that matches this contrast
                
                repcontrast = repmat(tmpindic(k,:),size(tmpindic2,1),1);    % replicate contrast
                if ~isempty(repcontrast)
                    [ib]=find(sum(abs(repcontrast - tmpindic2),2) == 0);
                    rmatx(e+k,i) = length(ib);
                else
                    % no contrasts even match any conditions for this study
                    % -- NaN?
                    rmatx(e+k,i) = 0;
                end
            end
        end
     
    end % loop through clusters
    
end % loop through studies

% make indicator for first fnames variable
%[indic,indinms] = string2indicator(fmatx);

%ustudyout = ustudyout';

ustudyout = DB.Study(allsuni);
    
return
