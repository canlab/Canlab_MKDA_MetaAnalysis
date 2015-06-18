function [Xi,alltasknms,condf,allti,testfield,prop_by_condition,num_by_condition] = Meta_Task_Indicator(DB,varargin)
% [Xi,alltasknms,condf,allti,testfield,prop_by_condition,num_by_condition] = Meta_Task_Indicator(DB,[data])
%
% Get task indicators and number / proportion of activating contrasts from DB
%
% needs to set up design:
%DB.(fields)   % lists of fields containing task conditions for each coordinate point
%DB.Contrast
%DB.connumbers
%
% Example:
% [Xi,alltasknms,condf,allti,testfield,prop_by_condition,num_by_condition]
% ... = Meta_Task_Indicator(DB,MC_Setup.unweighted_study_data);
%
% tor wager
% Created: 8/23/06

global dolist
dolist = 1;

if length(varargin) > 0,
    data = varargin{1};
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
        
        allti = [allti ti];
        alltasknms = [alltasknms tasknms];

    end
    
end
fprintf(1,'\n'); fprintf(1,'\n');

% -----------------------------------------------------
% allti is coords x tasks.  Convert to contrasts x tasks
% -----------------------------------------------------

ntasks = size(allti,2);
ncons = length(DB.connumbers);
Xi = zeros(ncons,ntasks);

for i = 1:ncons
    % which points in this contrast
    wh = DB.Contrast == DB.connumbers(i); 

    % indicators for this study
    ti = allti(wh,:);
    
    % resolve conflicts by taking the modal values
    % conflicts may be a sign of probs with database
    Xi(i,:) = mode(ti);

end

    
% return condition function (dummy codes) for chi-square test
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

% these must sum to 1 !  (they will be re-normed in Meta_Logistic, but
% summing to 1 is the standard weighting used in other parts of the
% software.)
DB.studyweight = w ./ sum(w);


% -----------------------------------------------------
% Get proportions by condition
% -----------------------------------------------------
prop_by_condition = [];
if exist('data','var')
    [icon,ctxtxi,betas,num_by_condition,prop_by_condition] = meta_apply_contrast(data,Xi,w,ones(1,size(Xi,2)));
end



return






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
    testfield = input('Enter field name or return if finished: ','s');

    if isempty(testfield)
        gook = 1;
    else
        gook = strmatch(testfield,N);
        if ~isempty(gook)
            gook = strcmp(testfield,N{gook}); 
            if ~gook, gook = []; end
        end %exact match
    end
    if isempty(gook), fprintf(1,'Not a field name.'); pause(1); fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b'); end
end

return





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

return

