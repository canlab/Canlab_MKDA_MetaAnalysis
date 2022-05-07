function DB = Meta_tables_of_variables(DB)
% Constructs some tables for each variable, saved in DB.TABLES, so that we can easily see the levels and numbers of points & contrasts for each variable
%
% DB = Meta_tables_of_variables(DB)
%
% - Must run Meta_Setup.m first
% - Assumes each study has one and only one unique text string assigned
% (no variation in capitalization, punctuation, etc. within a study)
%
% See also:
% [num, wh, weighted_num, conindx] = meta_count_contrasts(DB, testfield, fieldvalue)
% 

% Tor Wager, 5/7/2022

% initialize tables, in case we have selected levels/variables in DB,
% so that we want to re-start.

DB.TABLES = [];

% Check study var and fix if needed
if ~isfield(DB, 'study')
    if isfield(DB, 'Study')
        DB.study = DB.Study;
        
    else
        error('Enter DB.study or DB.Study');
    end
end

names = fieldnames(DB);

% exclude a few uninteresting variables
names(strcmp(names, 'x')) = [];
names(strcmp(names, 'y')) = [];
names(strcmp(names, 'z')) = [];
names(strcmp(names, 'xyz')) = [];

len = length(DB.x);  % number of points

for i = 1:length(names)
    
    fieldn = names{i};
    
    if size(DB.(fieldn), 1) == len % for variables with point annotations...
        
        DB.TABLES.(fieldn) = create_variable_table(DB, fieldn);
        
        % assume some vars not of interest for printing tables,
        % if n levels >= n contrasts
        
        if size(DB.TABLES.(fieldn), 1) < length(DB.pointind)
            
            disp(fieldn)
            disp(DB.TABLES.(fieldn))
            
            % check
            np = sum(DB.TABLES.(fieldn).NumPoints);
            nc = sum(DB.TABLES.(fieldn).NumContrasts);
            ns = sum(DB.TABLES.(fieldn).NumStudies);
            
            tt = table(np, nc, ns, 'VariableNames', {'TotalPoints' 'TotalContrasts' 'TotalStudies'});
            disp(tt);
            
            fprintf('Compare to: Total %3.0f points, %3.0f contrasts, and %3.0f studies in DB\n', length(DB.x), length(DB.pointind), length(unique(DB.study)));
            disp('(Studies may report multiple contrasts total may be > than num of unique studies.)')
            disp(' ')
            
        end % print
        
    end % add table
    
end % add all tables

disp('Added variable tables to DB.TABLES')

end % function





function t = create_variable_table(DB, fieldn)
% input: DB and field name (fieldn)
% output: table summarizing levels, points, and contrasts

if iscell(DB.(fieldn))
    % if cell array with text entries
    u = unique(DB.(fieldn));
    
else
    % assume numeric vector or matrix of column vectors
    u = unique(DB.(fieldn), 'rows');
    
end


clear npts ncons nstudies

for i = 1:size(u, 1)
    
    if iscell(DB.(fieldn))
        % if text
        wh = strcmp(DB.(fieldn), u{i});
    else
        % assume numeric vector or matrix of column vectors
        wh = ismember(DB.(fieldn), u(i, :), 'rows');
    end
    
    % count points
    npts(i, 1) = sum(wh);
    
    % count contrasts
    ncons(i, 1) = sum(wh(DB.pointind));
    
    
    % count studies
    study_list = DB.study(wh);

    nstudies(i, 1) = length(unique(study_list));
    
end

t = table(u, npts, ncons, nstudies, 'VariableNames', {'Level' 'NumPoints' 'NumContrasts' 'NumStudies'});

t.Properties.Description = fieldn;



end

