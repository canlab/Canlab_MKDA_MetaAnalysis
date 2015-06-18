function OUT = get_contrast_indicator(DB,testfield,varargin)
% OUT = get_contrast_indicator(DB,testfield,varargin)
% OUT = get_contrast_indicator(DB,'Method','create')
% OUT = get_contrast_indicator(DB,'valence','load',OUT)
% OUT = get_contrast_indicator(DB,'valence','create',[],{'pos' 'neg'})
% OUT = get_contrast_indicator(DB,'method','load',OUT,{'visual' 'auditory' 'recall' 'imagery'});
% Gets a structure with indicators and other stuff for density2 analysis
% Used in dbcontrast2density
%
% There are 2 modes this function runs in: 'load' and 'create'
% Specify this with a 3rd optional input argument
% Default is create, which creates a new set of contrast images and stores
% info in *info.mat files
%
% 'Load' looks for existing info.mat files and loads info about task
% conditions, etc. from those.  That way, you can use the same set of
% contrast images to conduct tesets on other variables, or compare
% specificity across different variables
% requires OUT as 2nd var argument; saved from dbcontrast2density
%
% This function returns sample sizes in indicator matrix
% -----------------------------------------------------
% identify independent contrasts
% -----------------------------------------------------

action = 'create';
if length(varargin) > 0, action = varargin{1};,end

switch action
    
    
    case 'load'
    %-Load function
    %----------------------------------------------------------------------- 
    OUT = varargin{2};
    [allconditions,pointcounts,studynames,allcondindic,allcondnames,conindex,PP] = loadindic(DB,OUT,testfield);    
        
    case 'create'
    %-Create function
    %-----------------------------------------------------------------------

% weight by sample size (opt input); requires Subjects field
% conindex contains indices for unique contrasts in vector of all peaks
% we need this to save variable info about each contrast
[allconditions,pointcounts,studynames,allcondindic,allcondnames,conindex] = dbcluster2indic(DB,DB,{testfield},1);


    otherwise
    %-Bad input
    %-----------------------------------------------------------------------
    error('Enter load or create')
end



fprintf(1,'Test variable is %s\n',testfield)
% number each independent contrast -- use to index contrast .img filenames
% throughout
connumbers = 1:length(studynames);  



%-----------------------------------------------------------------------
% Remove levels (columns) that were not specified, so we don't make maps
% for those
% This step does not remove contrasts from the analysis
% Just means we won't look for specificity for removed levels
% always remove NaN levels
%-----------------------------------------------------------------------

if length(varargin) > 2, 
    savelevels = varargin{3};,
    whsave = zeros(size(allcondnames));
else,   % save everything but NaNs
    savelevels = [];, 
    whsave = ones(size(allcondnames));
end

for i = 1:length(allcondnames), 
    if strcmp('NaN',allcondnames{i}), whomit(i)=1;, else,whomit(i)=0;,end, 

    if ~isempty(savelevels)
        for j = 1:length(savelevels)
            if strcmp(savelevels{j},allcondnames{i}), whsave(i)=1;,end, 
        end
    end
end
whsave(find(whomit)) = 0;
whsave = find(whsave);
allcondnames = allcondnames(whsave);
allcondindic = allcondindic(:,whsave);
% we would remove contrasts with none of the tasks here if we wanted to
% depends if you want to compute probs over all contrasts or only eligible
% ones
wh = find(all(allcondindic==0,2));
allconditions(wh) = [];
studynames(wh) = [];
pointcounts(wh,:) = [];
allcondindic(wh,:) = [];
connumbers(wh) = [];
conindex(wh) = [];
if exist('PP') == 1, PP(wh,:) = [];, OUT.PP = PP; end
fprintf(1,'%3.0f Contrasts are not eligible tasks: Removing these.\n',length(wh))
wh = [];

if size(studynames,1) == 1, studynames = studynames';, end

% this should now be redundant
% remove NaN contrasts
%wh = [];
%for i = 1:length(allconditions)
%    if isnan(allconditions{i}) | strcmp(allconditions{i},'NaN'),
%        wh(end+1) = i;
%    end
%end
%allconditions(wh) = [];
%studynames(wh) = [];
%pointcounts(wh,:) = [];
%allcondindic(wh,:) = [];
%connumbers(wh) = [];
%conindex(wh) = [];
%fprintf(1,'%3.0f Contrasts have NaN values for test field: Removing these.\n',length(wh))
%wh = [];
%for i = 1:length(allcondnames)
%    if isnan(allconditions{i}),
%        wh(end+1) = i;
%    end
%end
%allcondnames(wh) = [];


xyz = [DB.x DB.y DB.z];
eval(['allcond = DB.' testfield ';']); 
try, study = DB.Study;, catch, study = DB.study;, end

fprintf(1,'%3.0f Independent contrasts.\n',length(studynames))
fprintf(1,'\n');

OUT.allconditions = allconditions; 
OUT.pointcounts = pointcounts;
OUT.studynames = studynames;
OUT.allcondindic = allcondindic;
OUT.allcondnames = allcondnames;
OUT.testfield = testfield;
OUT.connumbers = connumbers;
OUT.conindex = conindex;


% -----------------------------------------------------
% Sample size weighting stuff (mostly output)
% -----------------------------------------------------

fprintf(1,'Weighting by sqrt of sample size is ON.\n');
fprintf(1,'Sample size statistics:\n');
fprintf(1,'\tTask frequencies\t');

task_count = sum(allcondindic > 0);           % number of contrasts for each task type
str = [repmat('%3.0f\t',1,size(allcondindic,2))];
fprintf(1,str,task_count)
fprintf(1,'\n')
avgn = nansum(allcondindic) ./ task_count;
fprintf(1,'\tAverage N\t');
str = [repmat('%3.2f\t',1,size(allcondindic,2))];
fprintf(1,str,avgn)
fprintf(1,'\n');
maxn = max(allcondindic);
fprintf(1,'\t Max N\t');
str = [repmat('%3.0f\t',1,size(allcondindic,2))];
fprintf(1,str,maxn)
fprintf(1,'\n');
minn = min(allcondindic(allcondindic>0));
fprintf(1,'\tMin N\t');
fprintf(1,str,minn)
fprintf(1,'\n');
fprintf(1,'\n');

% take square root for weighting
% get avg sqrt(n) to normalize scaling so relative sample sizes are used
allcondindic = allcondindic .^ .5;
avgn = nansum(allcondindic) ./ task_count;

for i = 1:length(studynames)
    n = max(allcondindic(i,:));     % sqrt(number of subjects in this contrast)
    wh_task = find(strcmp(allconditions{i},allcondnames));
    rootn(i) = n;
    try
        studyweight(i) = n ./ avgn(wh_task);  % sample size weighting by relative sample size
    catch
        studyweight(i) = NaN;   % if NaN allconditions
    end
end

% save sqrt(n) and weight value to apply to contrasts
% weight is sqrt(n) / avg n within task class
OUT.rootn = rootn;
OUT.studyweight = studyweight;


pause(3)

fprintf(1,'Contrasts and task levels (entries are sqrt of sample size)\n');
task_indicator_table(studynames,allconditions,allcondindic,allcondnames)

return









function [allconditions,pointcounts,studynames,allcondindic,allcondnames,conindex,PP] = loadindic(DB,OUT,testfield)

str = ['[levels,i3,j] = unique(DB.' testfield ');']; eval(str)

% preserve order
[i2,wh]=sort(i3); levels=levels(wh);

allcondnames = levels';

% -----------------------------------------------------
% get and check mat file names
% -----------------------------------------------------
disp('Loading *info.mat files: Loading file names');
% studynames = OUT.studynames;  don't do this, rather, re-build study names
% from ALL mat files; some studynames could be eliminated in some OUT
% structures.

    % d has ALL mat files, sorted
    d = dir(['*info.mat']);
    d = str2mat(d.name);
    d = sort_image_filenames(d);
    
for i = 1:size(d,1)
    
    % find name before first _ -- that's studyname
    %wh = find(d(i,:) == '_'); wh = wh(1);
    
    %studynames{i} = d(i:wh-1);
    
    %if findstr(studynames{i},d(i,:))
        % we're OK
        
    %else
    %    disp(['Warning!!! ' studynames{i} ' does not appear to be in correct order in file list.'])
    %end

end
    
% -----------------------------------------------------
% load mat files and return info
% -----------------------------------------------------
for i = 1:size(d,1)
    
    load(d(i,:))
    
    % find name before first _ -- that's studyname
    wh = find(d(i,:) == '_'); wh = wh(1);
    studynames{i} = d(i,1:wh-1);
    
    str = ['mycon = CONTRAST.' testfield '{1};']; eval(str);
    
    allconditions{i,1} = mycon;
    
    allcondindic(i,1:length(levels)) = 0;
    wh = find(strcmp(mycon,levels));
    allcondindic(i,wh) = CONTRAST.Subjects;
    
    % contrast number
    tmp = nums_from_text(d(i,:));
    tmp(tmp ~= real(tmp)) = [];
    conindex(i,1) = tmp(end);
    
    % point counts (we lose this info)
    pointcounts(i,1) = 1;
    
    % image file name 
    str = [studynames{i} '_contrast_' num2str(conindex(i)) '.img'];    % name of image
    if ~(exist(str) == 2), 
        fprintf(1,['Warning! Image %s does not exist.\n'],str);,
        
        % try to look for 2nd blank
            wh = find(d(i,:) == '_'); wh = wh(2);
            studynames{i} = d(i,1:wh-1);
            str = [studynames{i} '_contrast_' num2str(conindex(i)) '.img'];    % name of image
            disp(['Trying ' str ' instead.'])
            if ~(exist(str) == 2), 
                fprintf(1,['Warning! Image %s does not exist.\n'],str);,
            end
            
    end
    if ~(exist('PP')==1), PP = str;, else, PP = str2mat(PP,str);,end
    
end

return

