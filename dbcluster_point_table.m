function dbcluster_point_table(cl)
% function dbcluster_point_table(cl)
% looks for Study or study field in clusters
% determines length, and looks for other fields of the same length
% uses specified fields in a particular order, then other fields

N = fieldnames(cl(1));

% reorder - reverse order, bottom rows are top of list

wh = find(strcmp(N,'alltask')); if ~isempty(wh),tmp = N(wh);N(wh) = []; N = [tmp;N];,end
wh = find(strcmp(N,'Subjects')); if ~isempty(wh),tmp = N(wh);N(wh) = []; N = [tmp;N];,end
wh = find(strcmp(N,'t_Value')); if ~isempty(wh),tmp = N(wh);N(wh) = []; N = [tmp;N];,end
wh = find(strcmp(N,'Zscore')); if ~isempty(wh),tmp = N(wh);N(wh) = []; N = [tmp;N];,end
wh = find(strcmp(N,'Z')); if ~isempty(wh),tmp = N(wh);N(wh) = []; N = [tmp;N];,end
wh = find(strcmp(N,'distance')); if ~isempty(wh),tmp = N(wh);N(wh) = []; N = [tmp;N];,end
wh = find(strcmp(N,'z')); if ~isempty(wh),tmp = N(wh);N(wh) = []; N = [tmp;N];,end
wh = find(strcmp(N,'y')); if ~isempty(wh),tmp = N(wh);N(wh) = []; N = [tmp;N];,end
wh = find(strcmp(N,'x')); if ~isempty(wh),tmp = N(wh);N(wh) = []; N = [tmp;N];,end
wh = find(strcmp(N,'Study')); if ~isempty(wh),tmp = N(wh);N(wh) = []; N = [tmp;N];,end
wh = find(strcmp(N,'study')); if ~isempty(wh),tmp = N(wh);N(wh) = []; N = [tmp;N];,end


for i = 1:length(cl)
    
    % print coordinates
    fprintf(1,'Coordinate: %3.0f %3.0f %3.0f\n',cl(i).mm_center(1), cl(i).mm_center(2),cl(i).mm_center(3));
    
    if isfield(cl,'Study'),study = cl(i).Study; L = length(cl(i).Study); end
    if isfield(cl,'study'),study = cl(i).Study; L = length(cl(i).Study); end
    
    if ~(exist('L')==1), disp('No Study or study variable in clusters.'), return;  end
    
    % print header
    for k = 1:length(N)
        
        % figure out if this field is the right length
        eval(['x = cl(i).' N{k} ';']);
        go = 0;
        if iscell(x), if length(x) == L, go = 1; end %tmp = x{i}; end
        elseif size(x,1) == L, go = 1; % tmp = x(i,1);
        end
        
        % print, if it is
        if go
            print_f(N,k,length(N));
        end
    end
    fprintf(1,'\n')
    
    
    
    % print body rows
    
    for j = 1:L

        for k = 1:length(N)
            eval(['tmp = cl(i).' N{k} ';'])
            print_f(tmp,j,L);
        end
        fprintf(1,'\n')
    end
    fprintf(1,'\n')
           
end




function print_f(x,i,L,varargin)

go = 0;
if iscell(x), if length(x) == L, go = 1; tmp = x{i}; end
elseif size(x,1) == L, go = 1; tmp = x(i,1);
end

if go
    if ~isstr(tmp), tmp = sprintf('%3.2f',tmp);,end
    
    fprintf(1,'%s\t',tmp);
end

return




