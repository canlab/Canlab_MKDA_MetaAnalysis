function [pt,st] = contingency_table(varargin)
% function [pt,st] = contingency_table(varargin)
% 
% makes 2-way contingency tables for pairs of variables
% varargin arguments are variables
% vars must be column cell array vectors containing strings
%
% Examples:
% global study
% study = DB.Study;
% [pt, st] = contingency_table(DB.Method, DB.Stimuli);

for i = 1:length(varargin)

    % convert to string, then to cell array, if necessary
    if ~iscell(varargin{i})
    elseif ~ischar(varargin{i}{1})
        varargin{i} = num2str(cat(1,varargin{i}));
    end
    
    if ~iscell(varargin{i})
        varargin{i} = mat2cell(varargin{i},ones(size(varargin{i})),1);
    end
    
    % define categories (unique values) of a variable
    cats{i} = unique(varargin{i});
end


for i = 1:length(cats)-1
    
    % make a new table for each input variable
    
    % define unique values (categories) for each pair of variables
    cat1 = cats{i};
    cat2 = cats{i+1};
    
    for j = 1:length(cat1)
        for k = 1:length(cat2)
            
            % pointcount, studylist, studycount
            [pc,sl,sc] = countemstudies(varargin{i},cat1{j},varargin{i+1},cat2{k});
            
            % pointtable, studytable
            pt{i}(j,k) = pc;
            st{i}(j,k) = sc;
            
        end
    end
end

% for SINGLE input only
if length(cats) == 1
    cat1 = cats{1};
    for j = 1:length(cat1)
        % pointcount, studansylist, studycount
            [pc,sl,sc] = countemstudies(varargin{1},cat1{j});
            
            % pointtable, studytable
            pt{1}(j,1) = pc;
            st{1}(j,1) = sc;
    end
    
    % print table
    fprintf(1,'\nCounting peaks\n')
    fstr1 = repmat('%s',1,length(cats{1}));
    for j = 1:length(cats{1})
        fprintf(1,'%s\t', cats{1}{j})
        fprintf(1,'%4.0f\t',pt{1}(j,1))
        fprintf(1,'\n')
    end
    
    fprintf(1,'\nCounting studies\n')
    fstr1 = repmat('%s',1,length(cats{1}));
    for j = 1:length(cats{1})
        fprintf(1,'%s\t', cats{1}{j})
        fprintf(1,'%4.0f\t',st{1}(j,1))
        fprintf(1,'\n')
    end
    
else
    % for MULTIPLE inputs   
    % print output tables

    print_table(pt,cats)
    print_table(st,cats)

end

return





function print_table(pt,cats)

for i = 1:length(pt)
    
    fstr1 = repmat('%s',1,length(cats{i}));
    fstr2 = repmat('%s',1,length(cats{i+1}));
    
    % header row
    fprintf(1,'\n\t')
    for j = 1:length(cats{i+1})
        fprintf(1,'%s\t',cats{i+1}{j})
    end
    fprintf(1,'\n')
    
    % header col and lines of data
    for j = 1:length(cats{i})
        fprintf(1,'%s\t', cats{i}{j})
        
        for k = 1:length(cats{i+1})
            fprintf(1,'%4.0f\t',pt{i}(j,k))
        end
        
        fprintf(1,'\n')
    end
end   
return