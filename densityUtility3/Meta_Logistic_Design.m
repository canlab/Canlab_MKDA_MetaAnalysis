function [X,Xnms,DB,Xi,alltasknms,condf,testfield,conweights] = Meta_Logistic_Design(DB,varargin)
% [X,Xnms,DB,Xi,Xinms,condf,testfield,conweights] = Meta_Logistic_Design(DB,[control strings])
%
% Set up logistic regression design matrix from DB
%
% needs to set up design:
%DB.(fields)   % lists of fields containing task conditions for each coordinate point
%DB.pointind   % indices of which coord points are in which unique
%                contrast
%DB.x          % x coordinates, just to get # of points
%
% Optional 'control strings'
% {'empty'}   : to specify empty contrasts for each condition
%
% Outputs:
% X, design matrix
% Xnms, names
% DB, database, 
%   ... modified weights (w) to reflect new contrasts included
%   ... added fields to mark included points and contrasts in original list
% Xi, indicator matrix of which contrasts are in which conditions in the
% analysis
% Xinms, names if Xi
%
% tor wager
% Modified: 1/22/06

global dolist
dolist = 1;
doempty = 0;

if length(varargin) > 0,
    controlflags = varargin{1};
    for i = 1:length(controlflags)
        if strcmp(controlflags{i},'empty'), doempty = 1; end
    end
end

% -----------------------------------------------------
% Select conditions and build design matrix
% -----------------------------------------------------

testfield = 'xxx';
X = []; Xnms = {}; allti = []; alltasknms = [];

while ~isempty(testfield)
    
    % get field name and save for later display (image)
    testfield = getfield(DB); 
    
    if isempty(testfield), 
        % done
    else
            
        testfield1 = testfield;
        
        % check to see if each contrast has only one level
        meta_check_contrasts(DB,testfield)
        
        % get unique levels
        levels = getlevels(DB,testfield);

        % get indicators for these levels
        [ti,tasknms] = string2indicator(DB.(testfield),levels);

        % get contrast coded regressors (getcontrast subfunction is
        % replaced by meta_enter_contrasts)
        [con,conname,contype,conweights] = meta_enter_contrasts(ti,tasknms);

        % plot image of contrasts
        meta_plot_contrasts(DB,testfield,levels);
        title('Before database pruning');
        
        for i = 1:length(conname), conname{i} = [testfield '_' conname{i}]; end
        
        allti = [allti ti];
        alltasknms = [alltasknms tasknms];
        
        X = [X con];
        Xnms = [Xnms conname];
    end
    
end
fprintf(1,'\n'); fprintf(1,'\n');

% -----------------------------------------------------
% Prune database based on in-analysis points (contrasts)
% -----------------------------------------------------

% get only unique contrast entries
include = sum(allti,2);
%inc_con = include(DB.pointind);


%DB = prune(DB,include,inc_con);     % makes new pointind; this excludes
%point-wise

DB = Meta_Prune(DB,include);        % this excludes contrast-wise

meta_plot_contrasts(DB,testfield1,levels);
title('After pruning (valid contrasts)');
        
% -----------------------------------------------------
% Prune design matrix
% -----------------------------------------------------
Xpts = X;       % original

% get rid of excluded points based on incomplete/multi-level contrasts
X = X(find(DB.included_from_original),:);

% select unique contrasts
X = X(DB.pointind,:);

 
% prune indicator matrix
Xi = allti;
Xi = Xi(find(DB.included_from_original),:);

% select unique contrasts
Xi = Xi(DB.pointind,:);

% add empty contrasts if specified
if doempty
    if strcmp(contype,'contrast')
        X = add_empty_cons(X,Xnms);
    elseif strcmp(contype,'dummy')
        Xi = add_empty_cons(Xi,alltasknms);
    end
end

% image design matrix
 tor_fig(1,2); s = get(0,'ScreenSize'); set(gcf,'Position',[s(3).*.5 s(4).*.5 600 600]); 
    tmp = X;  set(gca,'YDir','Reverse');
    imagesc(tmp); colormap gray; title('Design matrix'); ylabel('Contrast'); xlabel('Regressors') 
    subplot(1,2,2); tmp = Xi; 
    set(gca,'YDir','Reverse');
    imagesc(tmp); colormap gray; title('Indicators'); ylabel('Contrast'); xlabel('Regressors') 
    drawnow
    
% end condition function (dummy codes) for chi-square test
[i,j] = find(Xi);
[i,s] = sort(i); 
condf = j(s);

% -----------------------------------------------------
% Re-do Final contrast weights
% -----------------------------------------------------

if isfield(DB,'SubjectiveWeights'),
    w = DB.rootn .* DB.SubjectiveWeights(DB.pointind);
else
    w = DB.rootn;
end

% if we have empty contrasts
% if empty images specified, ncontrasts will be greater than nimgs
if doempty
    addtow = zeros(size(Xi,1) - length(w),1);
    w = [w; addtow];
end

% these must sum to 1 !  (they will be re-normed in Meta_Logistic, but
% summing to 1 is the standard weighting used in other parts of the
% software.)
DB.studyweight = w ./ sum(w);



diary DESIGN_REPORT.txt

fprintf(1,'\nDesign Reporting\n---------------------------\nRegressors:\n')
for i = 1:length(Xnms), fprintf(1,'%3.0f\t%s\t%3.0f Contrasts (rows)\n',i,Xnms{i},sum(X(:,i)));end

empty = length(find(~sum(X,2)));
fprintf(1,'Intercept: %3.0f Contrast (sum of rows)\n',empty);

c = cond(X);
fprintf(1,'\nOrthogonality: \n\nCondition number (1 is good, higher is worse): %3.2f\n', c);

fprintf(1,'\nDependence (X''X) Matrix\n');
try
    print_matrix(X'*X,Xnms);
catch
end
fprintf(1,'\n')

diary off


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
        if len == length(DB.x), fprintf(1,'%s\n',N{i});end
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
    if isempty(gook), fprintf(1,'Not a field name.'); pause(1); fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b'); end
end

end





function levels = getlevels(DB,testfield);

eval(['levels = DB.' testfield ';']);
levels = levels(DB.pointind);
[pts,ind] = unique(levels);

fprintf(1,'Unique values of this variable:\n');

for i = 1:size(pts,1), 
    [num] = meta_count_contrasts(DB,testfield,pts{i});
    fprintf(1,'%3.0f\t%s\t%3.0f Contrasts\n',i,pts{i},num);
end
wh = input('Enter vector of levels to use: ');
levels = pts(wh);

end






%function [con,conname,type] = getcontrasts(ti,tasknms)


% type = 'contrast';
% 
% if size(ti,2) == 1
%     % simple predictor
%     con = ti(:,1);
%     conname = tasknms(1);
% elseif size(ti,2) == 2
%     % 2 levels, simple contrast
%     con = ti(:,1) - ti(:,2);
%     conname = {[tasknms{2} '-' tasknms{1}]};
% elseif size(ti,2) > 2
% 
%     if size(ti,2) == 3
%         type = input('Enter coding type: contrast vs. dummy:','s');
%     else
%         type = 'dummy';
%     end
% 
%     switch type
%         case 'contrast'
%             c = [1 -2 1; 1 1 -2]';
%             con = ti * c;
%             conname = {[tasknms{1} '+' tasknms{3} '-' tasknms{2}] [tasknms{1} '+' tasknms{2} '-' tasknms{3}]};
%         case 'dummy'
%             con = ti(:,1:end-1);
%             conname = tasknms(1:end-1);
% 
%             %wh = find(ti(:,end));  % unnecessary for dummy coding
%             %con(wh,:) = -1;
%             for i = 1:length(conname),conname{i} = [conname{i} '-' tasknms{end}];end
%         otherwise
%             error('type contrast or dummy!')
%     end
% else
%     % This should never happen.
%     %h = hadamard(4);
%     %error('only two or three levels allowed until program is expanded.');
%     % just treat this like a set of preds; don't re-parameterize
%     %con = ti(:,1:end-1);
%     %conname = tasknms(1:end-1);
%     %for i = 1:length(conname),conname{i} = [conname{i} '-' tasknms{end}];end
% end
% 


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
    if ismatrix(tmp), len = size(tmp,1);end
    
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
    if isempty(levels), levels = getlevels(DB,testfield); end

    % get indicators for these levels
    [ti,tasknms] = string2indicator(DB.(testfield),levels);


    tor_fig; s = get(0,'ScreenSize'); set(gcf,'Position',[s(3).*.5 s(4).*.5 600 600]); 
    tmp = [ti*100 DB.Contrast]; tmp(tmp == 0) = NaN;
    imagesc(tmp); set(gca,'XTick',1:size(tmp,2),'XTickLabel',[tasknms {'Contrast'}]); colormap prism
    cm = colormap(jet(size(tmp,1))); cm = cm(randperm(length(cm))',:);
    cm(1,:) = [1 1 1]; cm(2,:) = [1 0 0];
    colormap(cm); hold on; 
    for i = 1:length(DB.Contrast),wh1 = find(DB.Contrast==i); 
        if ~isempty(wh1)
            wh1 = wh1(end)+.5; plot([0 size(tmp,2)+.5], [wh1 wh1],'k'); 
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


function Xi = add_empty_cons(Xi,alltasknms)
[m,k]=size(Xi);
for i = 1:k
    n = input(['Enter number of empty contrasts for ' num2str(alltasknms{i}) ': ']);
    z = zeros(n,k); z(:,i) = 1;
    Xi = [Xi; z];
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


