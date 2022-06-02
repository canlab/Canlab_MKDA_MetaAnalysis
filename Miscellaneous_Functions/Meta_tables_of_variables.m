function DB = Meta_tables_of_variables(DB)
% Constructs some tables for each variable, saved in DB.TABLES, so that we can easily see the levels and numbers of points & contrasts for each variable
%
% DB = Meta_tables_of_variables(DB)
%
% - Must run Meta_Setup.m first
% - Assumes each study has one and only one unique text string assigned
% (no variation in capitalization, punctuation, etc. within a study)
%
% - Also extracts and adds Year field; if it doesn't exist yet, attempts to
% extract from .study
%
% - This function enforces some naming regularity not previously enforced,
% so is a good idea to run.  e.g., .study, .Year.
%
% Afterwards, you can do this kind of thing:
%     N_by_year_table = varfun(@mean, DB.TABLES.study_table, 'GroupingVariables','Year', 'InputVariables', 'mean_Subjects', 'OutputFormat','table');
%     ste_table = varfun(@ste, DB.TABLES.study_table, 'GroupingVariables','Year', 'InputVariables', 'mean_Subjects', 'OutputFormat','table');
%     
%     N_by_year_table = join(N_by_year_table, ste_table);
%     
%     create_figure('Sample size'); 
%     plot(DB.TABLES.study_table.Year, DB.TABLES.study_table.mean_Subjects, 'ko');
%     plot(N_by_year_table.Year, N_by_year_table.mean_mean_Subjects, 'LineWidth', 3, 'Color', 'b');
%     
%     xlabel('Year')
%     ylabel('Mean sample size per contrast')
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

%% Add Year

DB = add_Year_field_and_table(DB);

disp('Added variable tables to DB.TABLES')
disp('For final study table see DB.TABLES.study_table, and DB.TABLES.ccontrast_table')

end % function



% ----------------------------------------------------------
% ----------------------------------------------------------

%% Subfunctions

% ----------------------------------------------------------
% ----------------------------------------------------------



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



function DB = add_Year_field_and_table(DB)

DB = enforce_DB_Year_field(DB);
% Enforce that if we have a year field, it's DB.Year, and numeric

DB = enforce_DB_Study_field(DB);
% Enforce that if we have a study field, it's DB.study, lower-case s
    
[study, studyindx] = unique(DB.study, 'stable');

if ~isfield(DB, 'Year')

    disp('Extracting DB.Year from DB.study')

    % old way: regexp
    %     yr = zeros(size(DB.study));
    %     for i = 1:length(DB.study)
    %         yr(i, 1) = nums_from_text(DB.study{i});
    %     end

    % Assume we have 2-digit or 4-digit years, and want to include only these in your match:
    S = regexp(DB.study, '\d{2,4}', 'match', 'once', 'forceCellOutput');
    S(cellfun(@isempty, S)) = {'NaN'};
    yr = cellfun(@str2num, S);

    % now add 1900 to any years > 80 & < 99, assuming this is 2-digit entry
    wh = yr > 80 & yr < 99;
    yr(wh) = yr(wh) + 1900;

    % add 2000 to early numbers assuming 2-digit in 2000s
    wh = yr >= 0 & yr < 50; % good till 2050
    yr(wh) = yr(wh) + 2000;

    DB.Year = yr;

    %     yeartable = table(study, yr(wh), 'VariableNames', {'Study' 'Year'});
    %     disp(yeartable)

end


% Custom overall Contrast table

DB.TABLES.contrast_table = table(DB.study(DB.pointind), DB.Year(DB.pointind), DB.Subjects(DB.pointind), DB.Contrast(DB.pointind), 'VariableNames', {'study' 'Year' 'Subjects' 'Contrast'});

% Add sample size to study table: Mean subjects per contrast
mean_N = varfun(@mean, DB.TABLES.contrast_table, 'GroupingVariables','study', 'InputVariables', 'Subjects', 'OutputFormat','table');
mean_N.Properties.VariableNames{2} = 'Num_Contrasts';

% Custom overall study table
% note: rows are in different order from contrast table, so we must use
% join method

DB.TABLES.study_table = table(DB.study(studyindx), DB.Year(studyindx), 'VariableNames', {'study' 'Year'});

DB.TABLES.study_table = join(DB.TABLES.study_table, mean_N);

end


function DB = enforce_DB_Study_field(DB)
% Enforce that if we have a study field, it's DB.study
% Also add DB.studyindx

if ~isfield(DB, 'study') && isfield(DB, 'Study')

    DB.study = DB.Study;

end

[~, DB.studyindx] = unique(DB.study, 'stable');  % first entry in each study

end


function DB = enforce_DB_Year_field(DB)
% Enforce that if we have a year field, it's DB.Year, and numeric


if ~isfield(DB, 'Year') && isfield(DB, 'year')

    DB.Year = enforce_numeric_year(DB.year);

elseif isfield(DB, 'Year')

    DB.Year = enforce_numeric_year(DB.Year);

end

end

function year = enforce_numeric_year(year)
if isnumeric(year)
    % do nothing

else
    % convert
    year = cellfun(@str2num, year, 'UniformOutput', false);

    % if we have missing/other text values, need to further convert:
    if iscell(year)
        wh = ~cellfun(@isnumeric, year) | cellfun(@isempty, year);
        year(wh) = {NaN};
        year = [year{:}]';
    end

    if ~iscolumn(year), year = year'; end

end
end

