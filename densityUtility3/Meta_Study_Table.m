function Meta_Study_Table(DB,varargin)
% function Meta_Study_Table(DB,['study'])
%
% Prints text table of all independent contrasts for export
%
% looks for Study or study field in DB (also takes clusters, cl)
% determines length, and looks for other fields of the same length
% uses specified fields in a particular order, then other fields
%
% Similar to dbcluster point table, but prints one row per contrast

dostudy = 0;
for i = 1:length(varargin)
    if isstr(varargin{i})
        switch varargin{i}
           
            % functional commands
            case 'study', dostudy = 1;
         
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

N = fieldnames(DB(1));

% reorder field names - reverse order, bottom rows are top of list

wh = find(strcmp(N,'alltask')); if ~isempty(wh),tmp = N(wh);N(wh) = []; N = [tmp;N]; end
wh = find(strcmp(N,'Subjects')); if ~isempty(wh),tmp = N(wh);N(wh) = []; N = [tmp;N]; end
wh = find(strcmp(N,'t_Value')); if ~isempty(wh),tmp = N(wh);N(wh) = []; N = [tmp;N]; end
wh = find(strcmp(N,'Zscore')); if ~isempty(wh),tmp = N(wh);N(wh) = []; N = [tmp;N]; end
wh = find(strcmp(N,'Z')); if ~isempty(wh),tmp = N(wh);N(wh) = []; N = [tmp;N]; end
%wh = find(strcmp(N,'distance')); if ~isempty(wh),tmp = N(wh);N(wh) = []; N = [tmp;N]; end
%wh = find(strcmp(N,'z')); if ~isempty(wh),tmp = N(wh);N(wh) = []; N = [tmp;N]; end
%wh = find(strcmp(N,'y')); if ~isempty(wh),tmp = N(wh);N(wh) = []; N = [tmp;N]; end
%wh = find(strcmp(N,'x')); if ~isempty(wh),tmp = N(wh);N(wh) = []; N = [tmp;N]; end

% these variables are created in this script from study or Study and are
% the first two
N = ['year'; N];
N = ['Studynames'; N];


if ~isfield(DB,'pointind')
    disp('You need DB.pointind field.  Try Meta_Setup first')
    return; 
end
    
% DB can be clusters or single-element structure.  DO for each element
% (structure)
for i = 1:length(DB)
    
    % setup: get names and years with nice formatting
    % -----------------------------------------------
    wh = find(strcmp(N,'Study')); if ~isempty(wh), Studynames = DB(i).Study;  end %tmp = N(wh);N(wh) = []; N = [tmp;N]; end
    wh = find(strcmp(N,'study')); if ~isempty(wh), Studynames = DB(i).study;  end %tmp = N(wh);N(wh) = []; N = [tmp;N]; end
    L = length(Studynames);
    DB.Studynames = cell(L,1);
    DB.year = cell(L,1);
    
    % print coordinates for DB or cluster
    % -----------------------------------------------
    if isfield(DB,'mm_center')
        % this is a clusters struct, print center
        fprintf(1,'Coordinate: %3.0f %3.0f %3.0f\n',DB(i).mm_center(1), DB(i).mm_center(2),DB(i).mm_center(3));
    end
    
    if isfield(DB,'Study'),study = DB(i).Study; L = length(DB(i).Study); end
    if isfield(DB,'study'),study = DB(i).Study; L = length(DB(i).Study); end
    
    if ~(exist('L')==1), disp('No Study or study variable in DB.'), return; end
    
    % print header
    % -----------------------------------------------
    for k = 1:length(N)
        
        % figure out if this field is the right length
        eval(['x = DB(i).' N{k} ';']);
        go = 0;
        if iscell(x), if length(x) == L, go = 1; tmp = x{i}; end 
        elseif size(x,1) == L, go = 1; tmp = x(i,1);
        end
        
        % print, if it is
        if go
            print_f(N,k,length(N));
        end
    end
    fprintf(1,'\n')
    
    % get indices of rows: contrasts or studies
    % -----------------------------------------------
    if dostudy
        [tmp,DB.studyind] = unique(DB.Study);
        myfield = 'studyind';
    else
        myfield = 'pointind';
    end
    
    % print body rows
    % -----------------------------------------------
    conindx = DB.(myfield);
    if size(conindx,1) ~= 1, conindx = conindx'; end
    
    for j = conindx

        % get name and year with nice formatting
        [DB(i).Studynames{j},DB(i).year{j}] = sep_name_year(Studynames{j});
        
        for k = 1:length(N)
            eval(['tmp = DB(i).' N{k} ';'])
            print_f(tmp,j,L);
        end
        fprintf(1,'\n')
    end
    fprintf(1,'\n')
           
end




function print_f(x,i,L)

go = 0;
if iscell(x), if length(x) == L, go = 1; tmp = x{i}; end
elseif size(x,1) == L, go = 1; tmp = x(i,1);
end

if go
    
    if ~isstr(tmp), tmp = sprintf('%3.2f',tmp); end
    
    fprintf(1,'%s\t',tmp);
end

return




function [study,yr] = sep_name_year(studyname)
    
    yr = 'missing'; 
    study = studyname;
    [nums,whnums] = nums_from_text(studyname);
    if isnan(nums), return, end
  
    yr = nums(1);
    
    study = studyname;
    study(whnums) = [];
    
    if isempty(study) %% all numbers? assume not year...
        yr = 0;
    else
        % capitalize and format
        study(1) = upper(study(1));
    end
    
    if yr < 1900
        if yr > 80 && yr <= 99
            yr = yr + 1900;
        elseif yr < 20
            yr = yr + 2000;
        end
    end
    
    yr = sprintf('%04d',yr);
    

    return
    