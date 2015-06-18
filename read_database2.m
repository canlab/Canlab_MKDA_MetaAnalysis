function DB = read_database2(name, savename, varargin)
% Reads a CSV file into a meta-analysis database struct format using Matlab's importdata
%
% Usage:
%   -------------------------------------------------------------------------
%   [list outputs here] = function_name(name, savename, ['CalculateContrasts'])
%  
%   - runs Meta_Setup with a default 10 mm smoothing kernel
%   - Limited checking of mandatory fields - this could be improved
%   - Option to re-generate contrast numbers with 'CalculateContrasts'
%
%   Author and copyright information:
%   -------------------------------------------------------------------------
%       Copyright (C) 2013 Tor Wager
%   
%       This program is free software: you can redistribute it and/or modify
%       it under the terms of the GNU General Public License as published by
%       the Free Software Foundation, either version 3 of the License, or
%       (at your option) any later version.
%   
%       This program is distributed in the hope that it will be useful,
%       but WITHOUT ANY WARRANTY; without even the implied warranty of
%       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%       GNU General Public License for more details.
%   
%       You should have received a copy of the GNU General Public License
%       along with this program.  If not, see <http://www.gnu.org/licenses/>.
%  
%   Inputs:
%   -------------------------------------------------------------------------
%   name                    Name of CSV file, e.g., 'my_spreadsheet.csv'
%   'CalculateContrasts'    Optional keyword re-generates DB.Contrasts
%
%   Outputs:
%   -------------------------------------------------------------------------
%   DB          A meta-analysis database structure; also saves to disk
%  
%   Examples:
%   -------------------------------------------------------------------------
%  
%   %Go to dir for the analysis, beause files will be saved here:
%   cd('/Users/tor/Dropbox/Working_items/T_Neuroimmune_placebo_review/NN_Brain_body_neuroimaging_coordinates/Analysis_11_29_13')
%   name = 'Brain_body_autonomic.csv';
%   savename = 'DB_autonomic';
%   read_database2(name, savename);
%
%   See also:
%   read_database, for the main previous function to read in coordinates.
%   

dat = importdata(name);
dat.names = dat.textdata(1, :)';

%% Reformat into DB structure

ctr = 1;
istext = false(length(dat.names), 1);

for i = 1:length(dat.names)
    mytext = dat.textdata(2:end, i);
    
    istext(i) = any(~cellfun(@isempty, mytext)); % ~isempty(mytext{1});
    
    %if all(cellfun(@isempty, mytext))
    %    disp(['dat.names{' num2str(i) '} is empty.']);
    
    if istext(i)
        DB.(dat.names{i}) = mytext;
    else
        
        % catch in case this is empty col with neither data nor text
        % if such columns exist, will misalign though...
        if ctr > size(dat.data, 2)
            % empty field - use last existing
            DB.(dat.names{i}) = dat.data(:, ctr-1);
        else
            DB.(dat.names{i}) = dat.data(:, ctr);
            ctr = ctr + 1;
        end
    end
end

DB.xyz = [DB.x DB.y DB.z];

if isfield(DB, 'N') && iscell(DB.N)
    DB.N = cellfun(@str2num, DB.N);
    %DB.N = cat(1, DB.N{:});
end

if isfield(DB, 'Subjects') && iscell(DB.Subjects)
    DB.Subjects = cellfun(@str2num, DB.Subjects);
    %DB.Subjects = cat(1, DB.Subjects{:});
end

% Some function use lower case, some upper, unfortunately...
if isfield(DB, 'study') && ~isfield(DB, 'Study')
    DB.Study = DB.study;
end

%% Re/assign contrasts if requested

if any(strcmp(varargin, 'CalculateContrasts'))
    
    disp('RE-ASSIGNING CONTRAST NUMBERS');
    % Assign contrast numbers
    
    DB.Contrast = repmat({' '}, length(DB.x), 1);
    
    
    [DB.Contrast, first_peak_in_each_con] = get_contrast_indicator_improved(DB, dat.names(istext));
else
    disp('Using existing contrast numbers');
end

%% Check for required fields
fn = {'Contrast' 'Study' 'x' 'y' 'z'};
for f = fn
    if ~isfield(DB, f{1})
        error(['You must enter a field named ' f{1}]);
    end
end

% if Contrasts is text

if isfield(DB, 'Contrast') && iscell(DB.Contrast)
    contrasts = cellfun(@str2num, DB.Contrast, 'UniformOutput', false);
    contrasts = cat(1, contrasts{:});
    DB.Contrast_old = DB.Contrast;
    DB.Contrast = contrasts;
end


%% Check and replace - NaNs in sample size, etc.

if isfield(DB, 'N')
    DB.N(isnan(DB.N)) = nanmean(DB.N);
end

if isfield(DB, 'Subjects')
    DB.Subjects(isnan(DB.Subjects)) = nanmean(DB.Subjects);
end

%% Run Meta_Setup

DB = Meta_Setup(DB);
pause(3); 
close

%% Save

str = scn_get_datetime; str = str(1:end-6);

save(savename, 'DB');

end



