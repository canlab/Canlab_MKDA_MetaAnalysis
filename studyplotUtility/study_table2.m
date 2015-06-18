function out = study_table(study,varargin)
% function out = study_table(study,varargin)
%
% prints a table of studies; input the study variable, and input columns, in order
% columns with ones and zeros are converted to X's or blanks.
% (not implemented yet; everything must be cell arrays of strings now)
%
% study is assumed to be author name followed by 2-letter year code, e.g., 01 for 2001
%

% define vector of unique studies and initial stuff
[ustudy,wh] = unique(study);
[year,study] = getYear(study);
str = [repmat('%s\t',1,length(varargin)+2) '\n'];   % would be more elegant

% print header row
fprintf(1,'Author\tYear\t');
for i = 1:length(varargin),
    names{i} = inputname(i+1);, fprintf(1,'%s\t',names{i});, 
    
    % get unique
    %varargin{i} = varargin{i}(wh);
    
    if ~iscell(varargin{i}),
        varargin{i} = mat2cell(varargin{i},ones(length(varargin{i})),1);, 
        for j = 1:length(varargin{i}), varargin{i}{j} = num2str(varargin{i}{j});, end
    end 
    
    varargin{i}(strcmp(varargin{i},'nan')) = {'N/A'};
    
    % change 1's and 0's to X's, put in cell array
    % skip the 1's thing for now.  
    
end
fprintf(1,'\n')

out = cell(1);
% print rows
for i = 1:length(study)
    out{i} = sprintf('%s\t%s\t',study{i},year{i});
    for j = 1:length(varargin)
        out{i} = [out{i} sprintf('%s\t',varargin{j}{i})];
    end
    %out{i} = [out{i} sprintf('\n')];
end
        
out = unique(out);
for i = 1:length(out), fprintf(1,'%s\n',out{i});,end

out = str2mat(out);
        
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

