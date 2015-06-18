function DB = Meta_Prune(DB,include)
% DB = Meta_Prune(DB,include)
% 
% Prunes database given an indicator vector of points (peaks) to include
% Includes only contrasts for which ALL peaks are included!
% Thus, if any peaks are excluded, the whole contrast is excluded.
%
% Called in Meta_Select_Contrasts and Meta_Logistic_Design
%
% Used


% exclude any contrast with 1 or more points that are excluded.


exclude_points = include < 1;

% get contrast numbers to exclude
ex_connums = unique(DB.Contrast(exclude_points));


% get larger set of exclude points that excludes contrast-by-contrast
exclude_points = [];
for i = 1:length(ex_connums)
    wh = find(DB.Contrast == ex_connums(i));
    exclude_points = [exclude_points; wh];
end
    
% get indices in pointind to exclude
% exclude these indices in anything of length CONTRASTS
ex_conind = [];
for i = 1:length(ex_connums)
    wh = find(DB.connumbers == ex_connums(i));
    ex_conind = [ex_conind; wh];
end

% include indices
include_points = ones(size(include)); include_points(exclude_points) = 0;
include_cons = ones(size(DB.connumbers)); include_cons(ex_conind) = 0;

DB = prunedb(DB,include_points,include_cons);

DB.included_from_original = include_points;
DB.included_cons = include_cons;

fprintf(1,'Selected: %3.0f points, %3.0f contrasts, %3.0f studies.\n',length(DB.x),length(DB.connumbers),length(unique(DB.study)));


return





function DB = prunedb(DB,include,inc_con)
% Empty or field name from DB

N = fieldnames(DB);
whp = find(include);
whc = find(inc_con);

% special for DB.pointind
pind = zeros(size(DB.x)); pind(DB.pointind) = 1;
pind = pind(whp);

for i = 1:length(N)
    
    tmp = DB.(N{i}); len = length(tmp);
    if isnumeric(tmp), len = size(tmp,1); end
    
    if len == length(include), 
        DB.(N{i}) = tmp(whp);
    elseif len == length(inc_con)
       if isnumeric(tmp)
           DB.(N{i}) = tmp(whc,:);
       else
           DB.(N{i}) = tmp(whc);
       end
    end
    
end

% re-make pointind
for i = 1:length(DB.connumbers)
    wh = find(DB.Contrast == DB.connumbers(i));
    DB.pointind(i) = wh(1);
end


return