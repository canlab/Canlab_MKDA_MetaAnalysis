function [DB] = Meta_Select_Contrasts(DB)
% [DB] = Meta_Select_Contrasts(DB)
%
% Set up logistic regression design matrix from DB
%
% needs to set up design:
%DB.(fields)   % lists of fields containing task conditions for each coordinate point
%DB.pointind   % indices of which coord points are in which unique
%                contrast
%DB.x          % x coordinates, just to get # of points

global dolist
dolist = 1;

testfield = 'xxx';
X = []; Xnms = {};

while ~isempty(testfield)
    
    % get field name
    testfield = getfield(DB);
    
    if isempty(testfield), 
        % done
    else
            
        testfield1 = testfield;
        
        % check to see if each contrast has only one level
        meta_check_contrasts(DB,testfield)
        
        % get unique levels
        levels = getlevels(DB,testfield);

        % plot image of contrasts
        meta_plot_contrasts(DB,testfield,levels);
        title('Before database pruning');
        
        % get indicators for these levels
        [ti,tasknms] = string2indicator(DB.(testfield),levels);

        X = [X ti];
        Xnms = [Xnms tasknms];
    end
    
end

fprintf(1,'\n'); fprintf(1,'\n');


% -----------------------------------------------------
% Prune database based on in-analysis points (contrasts)
% -----------------------------------------------------

% get only unique contrast entries
include = sum(X,2);

DB = Meta_Prune(DB,include);        % this excludes contrast-wise

meta_plot_contrasts(DB,testfield1,levels);
title('After pruning (valid contrasts)');


% -----------------------------------------------------
% Re-do Final contrast weights
% -----------------------------------------------------

if isfield(DB,'SubjectiveWeights'),
    w = DB.rootn .* DB.SubjectiveWeights(DB.pointind);
else
    w = DB.rootn;
end

% these must sum to 1 !
DB.studyweight = w ./ sum(w);


end





function testfield = getfield(DB)
% Empty or field name from DB
global dolist

% list field names
fprintf(1,'Database field names\n---------------------------\n')
N = fieldnames(DB);
if dolist
    for i = 1:length(N),
        tmp = DB.(N{i}); len = length(tmp);
        if len == length(DB.x), fprintf(1,'%s\n',N{i});,end
    end

    fprintf(1,'\n');
    dolist = 0;
end

gook = [];
while isempty(gook)
    testfield = input('Enter field name or end if finished: ','s');

    if isempty(testfield),
        gook = 1;
    else
        gook = strmatch(testfield,N);
    end
    if isempty(gook), fprintf(1,'Not a field name.'); pause(1); fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');, end
end

end





function levels = getlevels(DB,testfield);

eval(['levels = DB.' testfield ';']);
levels = levels(DB.pointind);
[pts,ind] = unique(levels);

fprintf(1,'Unique values of this variable:\n');

for i = 1:size(pts,1), fprintf(1,'%3.0f\t%s\n',i,pts{i}); end
wh = input('Enter vector of levels to use: ');
levels = pts(wh);

end




function DB = prune(DB,include,inc_con)
% Empty or field name from DB

N = fieldnames(DB);
whp = find(include);
whc = find(inc_con);

% special for DB.pointind
%pind = zeros(size(DB.x)); pind(DB.pointind) = 1;
%pind = pind(whp);

for i = 1:length(N),
    
    tmp = DB.(N{i}); len = length(tmp);
    if ismatrix(tmp), len = size(tmp,1);,end
    
    if len == length(include), 
        DB.(N{i}) = tmp(whp);
    elseif len == length(inc_con), 
       if ismatrix(tmp),
           DB.(N{i}) = tmp(whc,:);
       else
           DB.(N{i}) = tmp(whc);
       end
    end
    
end

%DB.pointind = find(pind)';
%[DB.connumbers,DB.pointind] = unique(DB.Contrast);

% re-make pointind
for i = 1:length(DB.connumbers)
    wh = find(DB.Contrast == DB.connumbers(i));
    DB.pointind(i) = wh(1);
end

end




function meta_plot_contrasts(DB,testfield,levels)

    % get unique levels
    if isempty(levels), levels = getlevels(DB,testfield);, end

    % get indicators for these levels
    [ti,tasknms] = string2indicator(DB.(testfield),levels);


    tor_fig; s = get(0,'ScreenSize'); set(gcf,'Position',[s(3).*.5 s(4).*.5 600 600]); 
    tmp = [ti*100 DB.Contrast]; tmp(tmp == 0) = NaN;
    imagesc(tmp); set(gca,'XTick',1:size(tmp,2),'XTickLabel',[tasknms {'Contrast'}]); colormap prism
    cm = colormap(jet(size(tmp,1))); cm = cm(randperm(length(cm))',:);
    cm(1,:) = [1 1 1]; cm(2,:) = [1 0 0];
    colormap(cm); hold on; 
    for i = 1:length(DB.Contrast),wh1 = find(DB.Contrast==i);, 
        if ~isempty(wh1)
            wh1 = wh1(end)+.5; plot([0 size(tmp,2)+.5], [wh1 wh1],'k');, 
        end
    end
    set(gca,'YDir','Reverse'); ylabel('Peak Coordinate index');

    % print names
    for i = 1:size(ti,2)
        tmp = unique(DB.study(find(ti(:,i))))';
        fprintf(1,'\nVariable: %s  Level: %s  %3.0f studies\n',testfield,tasknms{i},length(tmp))
        fprintf(1,'%s\t',tmp{:})
        fprintf(1,'\n');
    end
        
end



% ISMATRIX: Returns 1 if the input matrix is 2+ dimensional, 0 if it is a scalar 
%           or vector.
%
%     Usage ismat = ismatrix(X)
%
% RE Strauss, 5/19/00

function ismat = ismatrix(X)
  [r,c] = size(X);
  if (r>1 && c>1)
    ismat = 1;
  else
    ismat = 0;
  end

end


