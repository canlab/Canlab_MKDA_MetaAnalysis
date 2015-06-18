function study_table(study,varargin)
% function study_table(study,varargin)
%
% prints a table of studies; input the study variable, and input columns, in order
% columns with ones and zeros are converted to X's or blanks.
% (not implemented yet; everything must be cell arrays of strings now)
%
% study is assumed to be author name followed by 2-letter year code, e.g., 01 for 2001
%

% define vector of unique studies and initial stuff
[ustudy,wh] = unique(study);
[year,ustudy] = getYear(ustudy);
str = [repmat('%s\t',1,length(varargin)+2) '\n'];   % would be more elegant

% print header row
fprintf(1,'Author\tYear\t');
for i = 1:length(varargin),
    names{i} = inputname(i+1);, fprintf(1,'%s\t',names{i});, 
    
    % get unique
    varargin{i} = varargin{i}(wh);
    
    if ~iscell(varargin{i}),
        varargin{i} = mat2cell(varargin{i},ones(length(varargin{i})),1);, 
        for j = 1:length(varargin{i}), varargin{i}{j} = num2str(varargin{i}{j});, end
    end 
    
    varargin{i}(strcmp(varargin{i},'nan')) = {'N/A'};
    
    % change 1's and 0's to X's, put in cell array
    % skip the 1's thing for now.  
    
end
fprintf(1,'\n')

% print rows
for i = 1:length(ustudy)
    fprintf(1,'%s\t%s\t',ustudy{i},year{i})
    for j = 1:length(varargin)
        fprintf(1,'%s\t',varargin{j}{i})
    end
    fprintf(1,'\n')
end
        
        
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



    % gender and stuff
    gen = unique(gender(strcmp(study,ustudy{i})));
    if any(strcmp(gen,'f')), c7 = 'X';, else c7 = ' ';,end
    if any(strcmp(gen,'m')), c8 = 'X';, else c8 = ' ';,end
    if any(strcmp(gen,'b')), c9 = 'X';, else c9 = ' ';,end
    
     val = unique(valence(strcmp(study,ustudy{i})));
    if any(strcmp(val,'pos')), c1 = 'X';, else c1 = ' ';,end
    if any(strcmp(val,'neg')), c2 = 'X';, else c2 = ' ';,end
    if any(strcmp(val,'mix')), c3 = 'X';, else c3 = ' ';,end
    
    app = unique(aav(strcmp(study,ustudy{i})));
    if any(strcmp(app,'approach')), c4 = 'X';, else c4 = ' ';,end
    if any(strcmp(app,'avoid')), c5 = 'X';, else c5 = ' ';,end
    if any(strcmp(app,'mix')), c6 = 'X';, else c6 = ' ';,end
    
    fprintf(1,'%s\t%s%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n',a{i},c,b{i}, ...
        c1,c2,c3,c4,c5,c6,c7,c8,c9)
end

